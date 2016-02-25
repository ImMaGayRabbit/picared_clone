    /*
     * The MIT License
     *
     * Copyright (c) 2015 The Broad Institute
     *
     * Permission is hereby granted, free of charge, to any person obtaining a copy
     * of this software and associated documentation files (the "Software"), to deal
     * in the Software without restriction, including without limitation the rights
     * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
     * copies of the Software, and to permit persons to whom the Software is
     * furnished to do so, subject to the following conditions:
     *
     * The above copyright notice and this permission notice shall be included in
     * all copies or substantial portions of the Software.
     *
     * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
     * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
     * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
     * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
     * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
     * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
     * THE SOFTWARE.
     */

    package picard.analysis;

    import htsjdk.samtools.SAMFileHeader;
    import htsjdk.samtools.SAMFileHeader.SortOrder;
    import htsjdk.samtools.SAMRecord;
    import htsjdk.samtools.SamReader;
    import htsjdk.samtools.SamReaderFactory;
    import htsjdk.samtools.reference.ReferenceSequence;
    import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
    import htsjdk.samtools.util.CloserUtil;
    import htsjdk.samtools.util.IOUtil;
    import htsjdk.samtools.util.Log;
    import htsjdk.samtools.util.ProgressLogger;
    import htsjdk.samtools.util.SequenceUtil;
    import picard.PicardException;
    import picard.cmdline.CommandLineProgram;
    import picard.cmdline.Option;
    import picard.cmdline.StandardOptionDefinitions;

    import java.io.File;
    import java.util.*;
    import java.util.concurrent.*;
    import java.util.concurrent.atomic.AtomicBoolean;
    import java.util.concurrent.locks.Lock;
    import java.util.concurrent.locks.ReentrantLock;

    /**
     * Super class that is designed to provide some consistent structure between subclasses that
     * simply iterate once over a coordinate sorted BAM and collect information from the records
     * as the go in order to produce some kind of output.
     *
     * @author Tim Fennell
     */
    public abstract class SinglePassSamProgram extends CommandLineProgram {
        @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
        public File INPUT;

        @Option(shortName = "O", doc = "File to write the output to.")
        public File OUTPUT;

        @Option(doc = "If true (default), then the sort order in the header file will be ignored.",
                shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
        public boolean ASSUME_SORTED = true;

        @Option(doc = "Stop after processing N reads, mainly for debugging.")
        public long STOP_AFTER = 0;

        private static final Log log = Log.getInstance(SinglePassSamProgram.class);

        /**
         * Final implementation of doWork() that checks and loads the input and optionally reference
         * sequence files and the runs the sublcass through the setup() acceptRead() and finish() steps.
         */
        @Override
        protected final int doWork() {
            makeItSo(INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, STOP_AFTER, Arrays.asList(this));
            return 0;
        }

        public static void makeItSo(final File input,
                                    final File referenceSequence,
                                    final boolean assumeSorted,
                                    final long stopAfter,
                                    final Collection<SinglePassSamProgram> programs) {

            // Setup the standard inputs
            IOUtil.assertFileIsReadable(input);
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(referenceSequence).open(input);

            // Optionally load up the reference sequence and double check sequence dictionaries
            final ReferenceSequenceFileWalker walker;
            if (referenceSequence == null) {
                walker = null;
            } else {
                IOUtil.assertFileIsReadable(referenceSequence);
                walker = new ReferenceSequenceFileWalker(referenceSequence);

                if (!in.getFileHeader().getSequenceDictionary().isEmpty()) {
                    SequenceUtil.assertSequenceDictionariesEqual(in.getFileHeader().getSequenceDictionary(),
                            walker.getSequenceDictionary());
                }
            }

            // Check on the sort order of the BAM file
            {
                final SortOrder sort = in.getFileHeader().getSortOrder();
                if (sort != SortOrder.coordinate) {
                    if (assumeSorted) {
                        log.warn("File reports sort order '" + sort + "', assuming it's coordinate sorted anyway.");
                    } else {
                        throw new PicardException("File " + input.getAbsolutePath() + " should be coordinate sorted but " +
                                "the header says the sort order is " + sort + ". If you believe the file " +
                                "to be coordinate sorted you may pass ASSUME_SORTED=true");
                    }
                }
            }

            // Call the abstract setup method!
            boolean anyUseNoRefReads = false;
            for (final SinglePassSamProgram program : programs) {
                program.setup(in.getFileHeader(), input);
                anyUseNoRefReads = anyUseNoRefReads || program.usesNoRefReads();
            }
            final ProgressLogger progress = new ProgressLogger(log);
            /*
            // Original Code:
            for (final SAMRecord rec : in) {
                // Setting Reference
                final ReferenceSequence ref;
                if (walker == null || rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                    ref = null;
                } else {
                    ref = walker.get(rec.getReferenceIndex());
                }
                //
                for (final SinglePassSamProgram program : programs) {
                    program.acceptRead(rec, ref);
                }

                progress.record(rec);

                // See if we need to terminate early?
                if (stopAfter > 0 && progress.getCount() >= stopAfter) {
                    break;
                }

                // And see if we're into the unmapped reads at the end
                if (!anyUseNoRefReads && rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                    break;
                }
            }
            */
            // --Setting up Executors--
            // Setting half of available processors to do the work
            // Maybe there is an optimal value for threads?
            // But we need at least 2 threads to do the job
            final ExecutorService executorService =
                    Executors.newFixedThreadPool(
                            Runtime.getRuntime().availableProcessors() < 4 ?
                                    2 :
                                    Runtime.getRuntime().availableProcessors()/2
                    );
            // --Variables
            final AtomicBoolean isStop = new AtomicBoolean(false);

            // --Constants--
            // Maybe there is some kind of algorithm to get ultimate list capacity like freeMemory/sizeOfChunk?
            final int LIST_CAPACITY = 10000;
            final int QUEUE_CAPACITY = 10;
//            final int SEM_CAPACITY = 10;

            // --Setting up some object stuff
            ArrayList<Object[]> pairs = new ArrayList<Object[]>(LIST_CAPACITY);
            final BlockingQueue<ArrayList<Object[]>> queue = new LinkedBlockingQueue<ArrayList<Object[]>>(QUEUE_CAPACITY);
            // No need for semaphore as queue restricting number of tasks running simultaneously
//            final Semaphore sem = new Semaphore(SEM_CAPACITY);

            final boolean finalAnyUseNoRefReads = anyUseNoRefReads;
            executorService.execute(new Runnable() {
                @Override
                public void run() {
                    while (true){
                        try {

                            final ArrayList<Object[]> pairsChunk = queue.take();
                            // Poison pill
                            if (pairsChunk.size() == 0){
                                return;
                            }
//                            sem.acquire();

                            executorService.execute(new Runnable() {
                                @Override
                                public void run() {
                                    for(final Object[] arr : pairsChunk){

                                        final SAMRecord rec = (SAMRecord) arr[0];
                                        final ReferenceSequence ref = (ReferenceSequence) arr[1];

                                        for (final SinglePassSamProgram program : programs) {
                                            program.acceptRead(rec, ref);
                                        }

                                        progress.record(rec);
                                        // See if we need to terminate early?
                                        if (stopAfter > 0 && progress.getCount() >= stopAfter) {
                                            isStop.set(true);
                                            return;
                                        }

                                        // And see if we're into the unmapped reads at the end
                                        if (!finalAnyUseNoRefReads && rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                                            isStop.set(true);
                                            return;
                                        }
                                    }
//                                    sem.release();
                                }
                            });
                        } catch (InterruptedException e) {
                            // Do nothing
                        }

                    }

                }
            });
            // Sometimes there are rare error, if exception is thrown and catched when iterating records
            // In that case ES isn't shut down, program won't exit. And that's kinda sad T__T
            for (final SAMRecord rec : in) {
                final ReferenceSequence ref;
                if (walker == null || rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                    ref = null;
                } else {
                    ref = walker.get(rec.getReferenceIndex());
                }
                // Checking if we need to stop
                if (isStop.get()){
                    // Shutting executorService
                    executorService.shutdownNow();
                    break;
                }
                pairs.add(new Object[]{rec, ref});
                if (pairs.size() < QUEUE_CAPACITY){
                    continue;
                }
                try {
                    queue.put(pairs);
                } catch (InterruptedException e) {
                    // Do nothing
                }
                pairs = new ArrayList<Object[]>(QUEUE_CAPACITY);
            }
            // This is not really good solution.
            // Have to redo this someday
            // when main thread is waiting for shutdown it won't bother
            // for a signal to stop all calculations
            // Maybe make some kind of Listener? Gotta think about it.
            if (!isStop.get()) {
                // Checking if array still has some pairs
                if (pairs.size() != 0) {
                    try {
                        queue.put(pairs);
                    } catch (InterruptedException e) {
                        // DO nothing
                    }
                }

                // Poison pill stuff
                pairs = new ArrayList<Object[]>();
                try {
                    queue.put(pairs);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }

                // Now shutting down Executor service
                executorService.shutdown();
            }
            // Just closing everything that is Closable
            CloserUtil.close(in);
            // There we can collect and compute final metrics (and write then in O-file? mb)
            for (final SinglePassSamProgram program : programs) {
                program.finish();
            }
        }

        /** Can be overriden and set to false if the section of unmapped reads at the end of the file isn't needed. */
        protected boolean usesNoRefReads() { return true; }

        /** Should be implemented by subclasses to do one-time initialization work. */
        protected abstract void setup(final SAMFileHeader header, final File samFile);

        /**
         * Should be implemented by subclasses to accept SAMRecords one at a time.
         * If the read has a reference sequence and a reference sequence file was supplied to the program
         * it will be passed as 'ref'. Otherwise 'ref' may be null.
         */
        protected abstract void acceptRead(final SAMRecord rec, final ReferenceSequence ref);

        /** Should be implemented by subclasses to do one-time finalization work. */
        protected abstract void finish();

    }
