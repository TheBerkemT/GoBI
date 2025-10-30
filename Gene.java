import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class Gene {
    //gene_id: given gene id on the .gtf file
    String gene_id;
    //gene_symbol: gene name given on the last part of the .gtf file
    String gene_symbol;
    String chr;
    String strand;
    // gene_transcripts: to link the transcripts to gene
    ArrayList<Transcript> gene_transcripts;

    public Gene(String geneId, String geneSymbol, String chromosome, String Gstrand){
        this.gene_id = geneId;
        this.gene_symbol = geneSymbol;
        this.chr = chromosome;
        this.strand = Gstrand;
        this.gene_transcripts = new ArrayList<Transcript>();
    }

    public void addTranscript(Transcript transcript){
        this.gene_transcripts.add(transcript);
    }

    public Transcript getTranscript(String transcriptId) {
        for (Transcript t : gene_transcripts) {
            if (t.transcript_id.equals(transcriptId)) return t;
        }
        return null;
    }

    public List<Transcript> getGene_transcripts() {
        return gene_transcripts;
    }

    //nports counter map for distinction
    public int countUniqueProteins() {
        HashSet<String> uniqueProteinIds = new HashSet<>();
        
        for (Transcript t : gene_transcripts) {
            for (String proteinId : t.protein_ids) {
                if (proteinId != null && !proteinId.isEmpty()) {
                    uniqueProteinIds.add(proteinId);
                }
            }
        }
        
        return uniqueProteinIds.size();
    }
}
