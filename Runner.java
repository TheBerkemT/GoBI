import java.io.*;
import java.util.*;

public class Runner {
    
    public static void main(String[] args) {
        System.out.println("NEW");
        // Check if no arguments provided - print usage
        if (args.length == 0) {
            printUsage();
            System.exit(0);
        }
        
        // Parse command line arguments
        String gtfFile = null;
        String outputFile = null;
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-gtf") && i + 1 < args.length) {
                gtfFile = args[i + 1];
                i++; // Skip next argument as it's the value
            } else if (args[i].equals("-o") && i + 1 < args.length) {
                outputFile = args[i + 1];
                i++; // Skip next argument as it's the value
            } else {
                System.err.println("Unknown parameter: " + args[i]);
                printUsage();
                System.exit(1);
            }
        }
        
        // Validate required parameters
        if (gtfFile == null || outputFile == null) {
            System.err.println("Error: Both -gtf and -o parameters are required.");
            printUsage();
            System.exit(1);
        }
        
        // Check if GTF file exists
        File gtfFileObj = new File(gtfFile);
        if (!gtfFileObj.exists()) {
            System.err.println("Error: GTF file not found: " + gtfFile);
            System.exit(1);
        }
        
        try {
            // Run the analysis
            System.out.println("=== EXON SKIPPING ANALYSIS ===");
            System.out.println("Input GTF: " + gtfFile);
            System.out.println("Output TSV: " + outputFile);
            System.out.println();
            
            long totalStartTime = System.currentTimeMillis();
            runAnalysis(gtfFile, outputFile);
            long totalEndTime = System.currentTimeMillis();
            
            System.out.println("\n=== ANALYSIS COMPLETE ===");
            System.out.println("Results written to: " + outputFile);
            System.out.println("Total execution time: " + formatTime(totalEndTime - totalStartTime));
            
        } catch (IOException e) {
            System.err.println("Error during analysis: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private static void runAnalysis(String gtfFile, String outputFile) throws IOException {
        // Step 1: Parse GTF file
        System.out.println("[Step 1/4] Parsing GTF file...");
        long startTime = System.currentTimeMillis();
        
        Map<String, Gene> genes = GTFParser.gtfToGeneMap(gtfFile);
        
        long parseTime = System.currentTimeMillis() - startTime;
        System.out.println("  ✓ Parsed " + genes.size() + " genes");
        System.out.println("  ✓ Time: " + formatTime(parseTime));
        System.out.println();
        
        // Step 2: Generate intron regions for all transcripts
        System.out.println("[Step 2/4] Generating intron regions...");
        startTime = System.currentTimeMillis();
        
        int totalTranscripts = 0;
        int totalIntrons = 0;
        int transcriptsWithIntrons = 0;
        
        for (Gene gene : genes.values()) {
            for (Transcript t : gene.gene_transcripts) {
                totalTranscripts++;
                ExonSkipNotice.intronRegioner(t);
                if (t.intron_regions.size() > 0) {
                    transcriptsWithIntrons++;
                }
                totalIntrons += t.intron_regions.size();
            }
        }
        
        long intronTime = System.currentTimeMillis() - startTime;
        System.out.println("  ✓ Processed " + totalTranscripts + " transcripts");
        System.out.println("  ✓ Transcripts with introns: " + transcriptsWithIntrons);
        System.out.println("  ✓ Generated " + totalIntrons + " introns");
        System.out.println("  ✓ Time: " + formatTime(intronTime));
        System.out.println();
        
        // Step 3: Detect exon skipping events
        System.out.println("[Step 3/4] Detecting exon skipping events...");
        startTime = System.currentTimeMillis();
        
        List<ExonSkipNotice> allEvents = new ArrayList<>();
        int genesWithEvents = 0;
        
        for (Gene gene : genes.values()) {
            List<ExonSkipNotice> geneEvents = ExonSkipNotice.exonSkips(gene);
            if (!geneEvents.isEmpty()) {
                genesWithEvents++;
                allEvents.addAll(geneEvents);
            }
        }
        
        long eventTime = System.currentTimeMillis() - startTime;
        System.out.println("  ✓ Found " + allEvents.size() + " exon skipping events");
        System.out.println("  ✓ Genes with events: " + genesWithEvents + " / " + genes.size());
        System.out.println("  ✓ Time: " + formatTime(eventTime));
        System.out.println();
        
        // Step 4: Write results to TSV
        System.out.println("[Step 4/4] Writing results to file...");
        startTime = System.currentTimeMillis();
        
        GTFParser.writeExonSkipsToTSV(allEvents, outputFile);
        
        long writeTime = System.currentTimeMillis() - startTime;
        System.out.println("  ✓ Output file: " + outputFile);
        System.out.println("  ✓ Time: " + formatTime(writeTime));
        
        // Summary
        System.out.println("\n--- TIMING SUMMARY ---");
        System.out.println("GTF Parsing:        " + formatTime(parseTime));
        System.out.println("Intron Generation:  " + formatTime(intronTime));
        System.out.println("Event Detection:    " + formatTime(eventTime));
        System.out.println("File Writing:       " + formatTime(writeTime));
        System.out.println("----------------------");
        System.out.println("Total:              " + formatTime(parseTime + intronTime + eventTime + writeTime));
    }
    
    private static String formatTime(long milliseconds) {
        if (milliseconds < 1000) {
            return milliseconds + " ms";
        } else if (milliseconds < 60000) {
            return String.format("%.2f s", milliseconds / 1000.0);
        } else {
            long minutes = milliseconds / 60000;
            long seconds = (milliseconds % 60000) / 1000;
            return String.format("%d min %d s", minutes, seconds);
        }
    }
    
    private static void printUsage() {
        System.out.println("USAGE:");
        System.out.println("  java -jar ExonSkipping.jar -gtf <GTF_file> -o <output_file>");
        System.out.println();
        System.out.println("REQUIRED PARAMETERS:");
        System.out.println("  -gtf <GTF_file>      Path to the input GTF annotation file");
        System.out.println("                       (Gene Transfer Format containing CDS features)");
        System.out.println();
        System.out.println("  -o <output_file>     Path to the output TSV file");
        System.out.println("                       (Tab-separated file with exon skipping events)");
        System.out.println();
        System.out.println("EXAMPLE:");
        System.out.println("  java -jar ExonSkipping.jar \\");
        System.out.println("       -gtf /path/to/gencode.v25.annotation.gtf \\");
        System.out.println("       -o results.tsv");
        System.out.println();
        System.out.println("DESCRIPTION:");
        System.out.println("  This tool analyzes GTF annotation files to identify exon skipping");
        System.out.println("  splicing events (ES-SE). It compares splice variants (SV) with");
        System.out.println("  wildtype (WT) transcripts to detect skipped exons.");
        System.out.println();
        System.out.println("OUTPUT FORMAT:");
        System.out.println("  The output TSV file contains the following columns:");
        System.out.println("  - Gene information (id, symbol, chromosome, strand)");
        System.out.println("  - Transcript counts (nprots, ntrans)");
        System.out.println("  - Intron coordinates (SV, WT)");
        System.out.println("  - Protein IDs (SV_prots, WT_prots)");
        System.out.println("  - Skipping statistics (min/max exons and bases)");
        System.out.println();
    }
}
