import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class GTFParser {

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