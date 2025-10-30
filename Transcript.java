import java.util.ArrayList;

public class Transcript {
    String transcript_id;
    ArrayList<String> protein_ids;
    ArrayList<CDS> exon_regions;
    
    //INTRON MUHABBET TRY 1
    ArrayList<IntronRegion> intron_regions;


    public void addIntron(IntronRegion intron){
        this.intron_regions.add(intron);
    }


    Transcript(String transcriptID){
        this.transcript_id = transcriptID;
        this.protein_ids = new ArrayList<>();
        this.exon_regions = new ArrayList<>();
        this.intron_regions = new ArrayList<>();
    }


}
