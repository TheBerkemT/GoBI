import java.util.ArrayList;
import java.util.Comparator;

import javax.swing.tree.TreeNode;

public class ExonSkipNotice {
    String geneID;
    String geneSymbol;   
    String chromosome;
    String strand;
    int nports;
    int ntrans;
    CDS SV_intron;              //the SV intron as start:end
    String WT_introns;          //all the WT introns within the SV intron separated by |as start:end
    String SV_prots;            //ids of the SV CDS, separated by |
    String WT_prots;            //ids of the WT CDS, separated by |
    int min_skipped_exon;       //the minimal number of skipped exons in any WT/SV pair
    int max_skipped_exon;       //the maximum number of skipped exons in any WT/SV pair
    int min_skipped_bases;      //the minimal number of skipped bases (joint length of skipped exons) in any WT/SV pair
    int max_skipped_bases;      //max skipped bases the maximum number of skipped bases (joint length of skipped exons) in any WT/SV pair.

    //We will derive the information from objekts made previously

    //ExonSkipNotice is our writer input objekt
    public ExonSkipNotice(String geneID, String geneSymbol, String chr, String strand, 
                          int nports, int ntrans, CDS SV_intron, String WT_introns, 
                          String SV_prots, String WT_prots, int min_skipped_exon, 
                          int max_skipped_exon, int min_skipped_bases, int max_skipped_bases){
        this.geneID = geneID;
        this.geneSymbol = geneSymbol;
        this.chromosome = chr;
    
    }

    public static void intronRegioner(Transcript t){
        //Benefit from ArrayList<> function .sort() 
        t.exon_regions.sort(Comparator.comparingInt(CDS::getStart));

        for(int k = 0; k < t.exon_regions.size()-1 ; k++){
            int intronStart = t.exon_regions.get(k).getEnd() + 1;
            int intronEnd = t.exon_regions.get(k + 1).getStart() - 1;
            CDS intron_region = new CDS(intronStart, intronEnd);
            
            t.addIntron(intron_region);
        }
    }

    /*
     * Determine the start and end regions of the intron region
     * will be matched on both the supposed Splice Variant and 
     * the Wild Type.
     * DELIVER AN ARRAYLIST AND IF ITS EMPTY THEN NO EXON SKIPS
     */
    public static boolean intronMatch(CDS svIntron, Transcript wt){
        boolean sameStart = false;
        boolean sameEnd = false;
        for(CDS intron : wt.intron_regions){
            if(svIntron.start == intron.start){
                sameStart = true;
            }
            if(svIntron.end == intron.end){
                sameEnd = true;
            }
        }
        return (sameStart && sameEnd);
    }
    
    //After intronMatch run this
    public static ArrayList<CDS> getRegionsIntrons(CDS sv_intron, Transcript t2){
        ArrayList<CDS> intronsInRegion = new ArrayList<>();

        for(CDS intron : t2.intron_regions){
            if(intron.start >= sv_intron.start && intron.end <= sv_intron.end){
                intronsInRegion.add(intron);
            }
        }
        return intronsInRegion;
    }
    /*
     *  IF THE ABOVE GIVEN intronsInRegion IS 1 
     *  THEN THE OTHER TRANSCRIPT IS ALSO A SV 
     *  AND SHOUDL BE ADDED TO SV_PROTS!!!!
     */
    

    //getting the skips CHANGE from arraylist to exonskipnotice OBJEKT
    //MIN MAX CALCULATIONS HERE MAYBE?
    //GET LENGTH HELPER FUNCTION?
    public static ArrayList<ExonSkipNotice> exonSkips(Gene geneToGetkips){
        ArrayList<ExonSkipNotice> exonSkipRegions = new ArrayList<>();

        for(Transcript sv : geneToGetkips.gene_transcripts){
            
            for(CDS sv_intron : sv.intron_regions){
                boolean wt_exists = false;
                //regions -> save for the skipping regions
                ArrayList<CDS> regions = new ArrayList<>();

                for(Transcript wt : geneToGetkips.gene_transcripts){
                    if(sv == wt) continue;
                
                    if(intronMatch(sv_intron, wt)){
                        regions = getRegionsIntrons(sv_intron, wt);
                        
                        //If regions.size()=1 then it is also a SV
                        //But there needs to be a WT to save them so check
                        if(regions.size() == 1){
                            //ADD TO THE SV_prots STRING!
                            
                        }
                        else if(regions.size() > 1){
                            wt_exists = true;
                            //ADD TO THE WT_prots STRING!
                            
                        }
                    }
                    if(wt_exists){
                        //IF WILD TYPE EXISTS THEN WE HAVE AT LEAST ONE SV-WT PAIR AND SHOULD SAVE IT

                        //EXONSKIPNOTICE OBJEKT HERE AND ADD IT TO THE ARRAYLIST<EXONSKIPOBJEKT>
                    }
                }
                
            }
        }
        return exonSkipRegions;
    }    

}
