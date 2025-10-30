import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class GTFParser {


    //TODO MOVE MAIN TO A RUNNER CLASS
    public static void main(String[] args) throws IOException {
        String gtfFile = "/Users/kazim.tamer/Desktop/Uni/WS25:26/GoBI/Exercise1/ServerGTF/gencode.v25.annotation.gtf";
        String outputFile = "test_output.tsv";
        
        System.out.println("=== EXON SKIPPING ANALYSIS TESTER ===\n");
        
        try {
            // Test 1: Parse GTF file
            System.out.println("TEST 1: Parsing GTF file...");
            long startTime = System.currentTimeMillis();
            Map<String, Gene> genes = GTFParser.gtfToGeneMap(gtfFile);
            long parseTime = System.currentTimeMillis() - startTime;
            
            System.out.println("✓ Successfully parsed GTF file");
            System.out.println("  - Total genes: " + genes.size());
            System.out.println("  - Parse time: " + parseTime + " ms");
            System.out.println();
            
            // Test 2: Check specific gene structure (SYNE1 - ENSG00000131018)
            System.out.println("TEST 2: Checking SYNE1 gene structure...");
            Gene syne1 = genes.get("ENSG00000131018");
            if (syne1 != null) {
                System.out.println("✓ Found SYNE1 gene");
                System.out.println("  - Gene ID: " + syne1.gene_id);
                System.out.println("  - Gene Symbol: " + syne1.gene_symbol);
                System.out.println("  - Chromosome: " + syne1.chr);
                System.out.println("  - Strand: " + syne1.strand);
                System.out.println("  - Number of transcripts: " + syne1.gene_transcripts.size());
                System.out.println("  - Number of unique proteins: " + syne1.countUniqueProteins());
                
                // Show first few transcripts
                System.out.println("  - Sample transcripts:");
                int count = 0;
                for (Transcript t : syne1.gene_transcripts) {
                    System.out.println("    * " + t.transcript_id + " (exons: " + t.exon_regions.size() + 
                                    ", proteins: " + t.protein_ids.size() + ")");
                    if (++count >= 3) break;
                }
            } else {
                System.out.println("✗ SYNE1 gene not found! Check gene_id extraction.");
            }
            System.out.println();
            
            // Test 3: Generate introns for all transcripts
            System.out.println("TEST 3: Generating intron regions...");
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
                        totalIntrons += t.intron_regions.size();
                    }
                }
            }
            long intronTime = System.currentTimeMillis() - startTime;
            System.out.println("✓ Generated introns");
            System.out.println("  - Total transcripts: " + totalTranscripts);
            System.out.println("  - Transcripts with introns: " + transcriptsWithIntrons);
            System.out.println("  - Total introns generated: " + totalIntrons);
            System.out.println("  - Time: " + intronTime + " ms");
            System.out.println();

            // Test 4: Exon skipping analysis
            System.out.println("TEST 4: Detecting exon skipping events...");
            startTime = System.currentTimeMillis();
            List<ExonSkipNotice> allEvents = new ArrayList<>();

            for (Gene gene : genes.values()) {
                List<ExonSkipNotice> geneEvents = ExonSkipNotice.exonSkips(gene);
                if (!geneEvents.isEmpty()) {
                    allEvents.addAll(geneEvents);
                }
            }

            long eventTime = System.currentTimeMillis() - startTime;
            System.out.println("✓ Exon skipping detection complete");
            System.out.println("  - Total exon skipping events: " + allEvents.size());
            System.out.println("  - Time: " + eventTime + " ms");
            System.out.println();

            // Test 5: Write results to TSV
            System.out.println("TEST 5: Writing results to TSV...");
            startTime = System.currentTimeMillis();
            GTFParser.writeExonSkipsToTSV(allEvents, outputFile);
            long writeTime = System.currentTimeMillis() - startTime;
            System.out.println("✓ TSV file written successfully");
            System.out.println("  - Output: " + outputFile);
            System.out.println("  - Time: " + writeTime + " ms");
            System.out.println();

            System.out.println("=== ANALYSIS COMPLETE ===");
            
        } catch (IOException e) {
            System.err.println("Error while reading GTF file: " + e.getMessage());
            e.printStackTrace();
        }
    
    }

    public static Map<String, Gene> gtfToGeneMap(String fileName) throws IOException {
        Map<String, Gene> genes = new HashMap<>();
    
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue; // Skip comment lines
            
                String[] columns = line.split("\t");
                String feature = columns[2];

                if (!feature.equals("CDS") && !feature.equals("transcript")) continue;

                String attributes = columns[8]; // The last column with gene_id, transcript_id, etc.
            
                // Extract the gene_id and transcript_id from the attributes string
                String geneId = extractAttribute(attributes, "gene_id");
                String transcriptId = extractAttribute(attributes, "transcript_id");
                String geneName = extractAttribute(attributes, "gene_name");
                String exon_id = extractAttribute(attributes, "exon_id");
                String protein_id = extractAttribute(attributes, "protein_id");
                
                String chromosome = columns[0];
                String strand = columns[6];

                // Get or create the Gene
                Gene gene = genes.get(geneId);
                if (gene == null) {
                    gene = new Gene(geneId, geneName, chromosome, strand);
                    genes.put(geneId, gene);
                }

                // Get or create the Transcript
                Transcript transcript = gene.getTranscript(transcriptId);
                if (transcript == null) {
                    transcript = new Transcript(transcriptId);
                    gene.addTranscript(transcript);
                }

                // Add CDS region
                CDS exon = new CDS(Integer.parseInt(columns[3]), Integer.parseInt(columns[4]),exon_id);
                transcript.exon_regions.add(exon);
                if (protein_id != null && !transcript.protein_ids.contains(protein_id))
                    transcript.protein_ids.add(protein_id);
            }
        }
        return genes;
    }


    private static String extractAttribute(String attributes, String key) {
        String[] keyValuePairs = attributes.split(";");
        for (String pair : keyValuePairs) {
            pair = pair.trim();
            if (pair.startsWith(key)) {
                return pair.split("\"")[1];
            }
        }
        return null;
    }


    public static void writeExonSkipsToTSV(List<ExonSkipNotice> events, String outFile) throws IOException {
        try (BufferedWriter w = new BufferedWriter(new FileWriter(outFile))) {
            w.write("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_prots\tSV_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases\n");
            for (ExonSkipNotice e : events) {

                String wtIntronsStr = e.WT_introns.stream()
                                                  .flatMap(List::stream)
                                                  .map(IntronRegion::toString)
                                                  .collect(Collectors.joining("|"));
                
                
                w.write(String.join("\t",
                    e.geneID,
                    e.geneSymbol,
                    e.chromosome,
                    e.strand,
                    String.valueOf(e.nports),
                    String.valueOf(e.ntrans),
                    e.SV_intron.toString(),
                    wtIntronsStr,
                    String.join("|", e.SV_prots),
                    String.join("|", e.WT_prots),
                    String.valueOf(e.min_skipped_exon),
                    String.valueOf(e.max_skipped_exon),
                    String.valueOf(e.min_skipped_bases),
                    String.valueOf(e.max_skipped_bases)
                ));
                w.newLine();
            }
        }
    }
}