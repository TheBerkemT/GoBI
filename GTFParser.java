import java.io.*;
import java.util.*;

public class GTFParser {


    //TODO MOVE MAIN TO A RUNNER CLASS
    public static void main(String[] args) throws IOException {
        String fileName = "/Users/kazim.tamer/Desktop/Uni/WS25:26/GoBI/Exercise1/ClOutputs/ENSG00000131018.input";
        Map<String, Gene> genes = gtfToGeneMap(fileName);
        

        // Print out genes and their transcripts
        for (Gene gene : genes.values()) {
            System.out.println(genes.size());
            System.out.println(genes.get(gene.gene_id).getGene_transcripts().size());
            System.out.println("Gene ID: " + gene.gene_id + ", Gene Name: " + gene.gene_symbol);
            for (Transcript transcript : gene.getGene_transcripts()) {
                System.out.println("Transcript ID: " + transcript.transcript_id);
                for(int i = 0; i < transcript.exon_regions.size(); i++){
                    System.out.println(i);
                    System.out.println(transcript.exon_regions.get(i).start);
                    System.out.println(transcript.exon_regions.get(i).end);
                }
                
            }
        }
        
        for(Gene gene : genes.values()){
            for(Transcript t : gene.gene_transcripts){
                ExonSkipNotice.intronRegioner(t);
                System.out.println(t.transcript_id);
            }
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
            
                // Get or create the Gene
                Gene gene = genes.get(geneId);
                if (gene == null) {
                    gene = new Gene(geneId, geneName);
                    genes.put(geneId, gene);
                }

                // Get or create the Transcript
                Transcript transcript = gene.getTranscript(transcriptId);
                if (transcript == null) {
                    transcript = new Transcript(transcriptId, geneId);
                    gene.addTranscript(transcript);
                }

                // Add CDS region
                CDS exon = new CDS(Integer.parseInt(columns[3]), Integer.parseInt(columns[4]));
                transcript.exon_regions.add(exon);
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
}