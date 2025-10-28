import java.util.ArrayList;

public class Transcript {
    String transcript_id;
    String protein_id;
    ArrayList<CDS> exon_regions;
    
    //INTRON MUHABBET TRY 1
    ArrayList<CDS> intron_regions;
    public void addIntron(CDS intron){
        this.intron_regions.add(intron);
    }


    Transcript(String transcriptID, String protID){
        this.transcript_id = transcriptID;
        this.protein_id = protID;
        this.exon_regions = new ArrayList<>();
        this.intron_regions = new ArrayList<>();
    }


}
