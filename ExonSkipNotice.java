import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;


public class ExonSkipNotice {
    String geneID;
    String geneSymbol;   
    String chromosome;
    String strand;
    int nports;
    int ntrans;
    IntronRegion SV_intron;              //the SV intron as start:end
    ArrayList<ArrayList<IntronRegion>> WT_introns;          //all the WT introns within the SV intron separated by |as start:end
    ArrayList<String> SV_prots;            //ids of the SV CDS, separated by |
    ArrayList<String> WT_prots;            //ids of the WT CDS, separated by |
    int min_skipped_exon;       //the minimal number of skipped exons in any WT/SV pair
    int max_skipped_exon;       //the maximum number of skipped exons in any WT/SV pair
    int min_skipped_bases;      //the minimal number of skipped bases (joint length of skipped exons) in any WT/SV pair
    int max_skipped_bases;      //max skipped bases the maximum number of skipped bases (joint length of skipped exons) in any WT/SV pair.


    //ExonSkipNotice is our writer input objekt                 Set used to be ArrayList on both constructor and
    public ExonSkipNotice(String geneID, String geneSymbol, String chr, String strand, 
                          int nports, int ntrans, IntronRegion SV_intron, ArrayList<ArrayList<IntronRegion>> WT_introns, 
                          ArrayList<String> WT_prots, ArrayList<String> SV_prots, int min_skipped_exon, 
                          int max_skipped_exon, int min_skipped_bases, int max_skipped_bases){
        this.geneID = geneID;
        this.geneSymbol = geneSymbol;
        this.chromosome = chr;
        this.strand = strand;
        this.nports = nports;
        this.ntrans = ntrans;
        this.SV_intron = SV_intron;
        this.WT_introns = WT_introns;
        this.SV_prots = SV_prots;
        this.WT_prots = WT_prots;
        this.min_skipped_exon = min_skipped_exon;
        this.max_skipped_exon = max_skipped_exon;
        this.min_skipped_bases = min_skipped_bases;
        this.max_skipped_bases = max_skipped_bases;
        
        }

    public static void intronRegioner(Transcript t){
        //Benefit from ArrayList<> function .sort() 
        t.exon_regions.sort(Comparator.comparingInt(CDS::getStart));

        for(int k = 0; k < t.exon_regions.size()-1 ; k++){
            int intronStart = t.exon_regions.get(k).getEnd() + 1;
            int intronEnd = t.exon_regions.get(k + 1).getStart();
            
            if (intronStart <= intronEnd) {
                IntronRegion intron_region = new IntronRegion(intronStart, intronEnd);
                t.addIntron(intron_region);
            } 
        }
    }

    /*
     * Determine the start and end regions of the intron region
     * will be matched on both the supposed Splice Variant and 
     * the Wild Type.
     * DELIVER AN ARRAYLIST AND IF ITS EMPTY THEN NO EXON SKIPS
     */
    public static boolean intronMatch(IntronRegion svIntron, Transcript wt){
        boolean sameStart = false;
        boolean sameEnd = false;
        for(IntronRegion intron : wt.intron_regions){
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
    public static ArrayList<IntronRegion> getRegionsIntrons(IntronRegion sv_intron, Transcript t2){
        ArrayList<IntronRegion> intronsInRegion = new ArrayList<>();

        for(IntronRegion intron : t2.intron_regions){
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
    
    public static ArrayList<ExonSkipNotice> exonSkips(Gene geneToGetSkips) {
        ArrayList<ExonSkipNotice> exonSkipRegions = new ArrayList<>();
        HashMap<String, ExonSkipNotice> skipMap = new HashMap<>();

        for (Transcript sv : geneToGetSkips.gene_transcripts) {
            for (IntronRegion sv_intron : sv.intron_regions) {

                boolean wt_exists = false;
                ArrayList<String> sv_temp = new ArrayList<>();
                sv_temp.add(protIdParser(sv.protein_ids));

                ArrayList<String> wt_temp = new ArrayList<>();
                ArrayList<ArrayList<IntronRegion>> wt_regions = new ArrayList<>();
                ArrayList<Transcript> wt_transcripts_list = new ArrayList<>();

                for (Transcript wt : geneToGetSkips.gene_transcripts) {
                    if (sv == wt) continue;

                    if (intronMatch(sv_intron, wt)) {
                        ArrayList<IntronRegion> regions = getRegionsIntrons(sv_intron, wt);

                        if (regions.size() == 1) {
                            sv_temp.add(protIdParser(wt.protein_ids));
                        }
                        else if (regions.size() > 1) {
                            wt_exists = true;
                            wt_regions.add(regions);
                            wt_transcripts_list.add(wt);
                            wt_temp.add(protIdParser(wt.protein_ids));
                        }
                    }
                }

                if (wt_exists) {
                    // DEDUPLICATE WT REGIONS
                    Set<String> seenPatterns = new HashSet<>();
                    ArrayList<ArrayList<IntronRegion>> uniqueWtRegions = new ArrayList<>();
                    ArrayList<Transcript> uniqueWtTranscripts = new ArrayList<>();
                    
                    for (int i = 0; i < wt_regions.size(); i++) {
                        ArrayList<IntronRegion> regionList = wt_regions.get(i);
                        
                        String pattern = regionList.stream()
                                                .map(IntronRegion::toString)
                                                .collect(Collectors.joining("|"));
                        
                        if (!seenPatterns.contains(pattern)) {
                            seenPatterns.add(pattern);
                            uniqueWtRegions.add(regionList);
                            uniqueWtTranscripts.add(wt_transcripts_list.get(i));
                        }
                    }
                    
                    // Calculate with deduplicated data
                    int[] exonmm = min_max_exons(uniqueWtRegions);
                    int[] basemm = min_max_bases(uniqueWtRegions, sv_intron, uniqueWtTranscripts);

                    int nprots = geneToGetSkips.countUniqueProteins();
                    int ntrans = geneToGetSkips.gene_transcripts.size();

                    String key = geneToGetSkips.gene_id + "_" + sv_intron.toString();
                    if (!skipMap.containsKey(key)) {
                        ExonSkipNotice event = new ExonSkipNotice(
                            geneToGetSkips.gene_id,
                            geneToGetSkips.gene_symbol,
                            geneToGetSkips.chr,
                            geneToGetSkips.strand,
                            nprots,
                            ntrans,
                            sv_intron,
                            uniqueWtRegions,  // Use deduplicated
                            sv_temp,          // SV_prots (will be swapped in constructor)
                            wt_temp,          // WT_prots (will be swapped in constructor)
                            exonmm[0],
                            exonmm[1],
                            basemm[0],
                            basemm[1]
                        );
                        skipMap.put(key, event);
                    }
                }
            }
        }

        exonSkipRegions.addAll(skipMap.values());
        return exonSkipRegions;
    }




    //YOU NEED A EXONSKIPNOTICE WRITER, IN IT DO START+:+END FOR INTRON REGIONS
    public static int[] min_max_bases(ArrayList<ArrayList<IntronRegion>> wt_regions, 
                                   IntronRegion sv_intron,
                                   ArrayList<Transcript> wt_transcripts) {
        int max = 0;
        int min = Integer.MAX_VALUE;

        // For each WT transcript that has this skipping pattern
        for (int i = 0; i < wt_transcripts.size(); i++) {
            Transcript wt = wt_transcripts.get(i);
            int totalBases = 0;
            
            // Find all exons in this WT that fall within the SV intron region
            for (CDS exon : wt.exon_regions) {
                // Check if exon is completely within the SV intron boundaries
                if (exon.start >= sv_intron.start && exon.end <= sv_intron.end) {
                    totalBases += (exon.end - exon.start + 1);
                }
                // Handle partial overlaps if needed (depends on your definition)
                // For now, we only count fully contained exons
            }

            if (totalBases > max) max = totalBases;
            if (totalBases < min) min = totalBases;
        }

        if (wt_transcripts.isEmpty()) {
            min = 0;
            max = 0;
        }

        return new int[]{min, max};
    }
    
    public static int[] min_max_exons(ArrayList<ArrayList<IntronRegion>> wt_regions) {
        int max = 0;
        int min = Integer.MAX_VALUE;

        for (ArrayList<IntronRegion> irl : wt_regions) {
            int size = irl.size();
            if (size > max) max = size;
            if (size < min) min = size;
        }

        // handle edge case: empty input
        if (wt_regions.isEmpty()) {
            min = 0;
            max = 0;
        } else {
            //adjust since we are counting introns, just do -1 since intron-exon-intron for skips
            if (min > 0) min = min - 1;
            if (max > 0) max = max - 1;
        }

        return new int[]{min, max};
    }

    //use this when writing into the tsv file
    public static String protIdParser(ArrayList<String> protIds){
        String s = "";
        for(String id : protIds){
            if(id==null) continue;
            
            if(s.length()==0){
                s = s + id;
            } else {
                s = s + "|" + id;
            }
        }
        return s;
    }
}


/*
 * public static ArrayList<ExonSkipNotice> exonSkips(Gene geneToGetSkips) {
        ArrayList<ExonSkipNotice> exonSkipRegions = new ArrayList<>();
        HashMap<String, ExonSkipNotice> skipMap = new HashMap<>();

        for (Transcript sv : geneToGetSkips.gene_transcripts) {
            // Loop over all SV introns
            for (IntronRegion sv_intron : sv.intron_regions) {

                boolean wt_exists = false;
                ArrayList<String> sv_temp = new ArrayList<>();
                sv_temp.add(protIdParser(sv.protein_ids));

                ArrayList<String> wt_temp = new ArrayList<>();
                ArrayList<ArrayList<IntronRegion>> wt_regions = new ArrayList<>();
                ArrayList<Transcript> wt_transcripts_list = new ArrayList<>();

                // Compare to all other transcripts in this gene
                for (Transcript wt : geneToGetSkips.gene_transcripts) {
                    if (sv == wt) continue;

                    // If WT transcript shares start/end boundaries
                    if (intronMatch(sv_intron, wt)) {
                        ArrayList<IntronRegion> regions = getRegionsIntrons(sv_intron, wt);

                        // Same intron → also SV
                        if (regions.size() == 1) {
                            sv_temp.add(protIdParser(wt.protein_ids));
                        }
                        // More than 1 intron inside → WT form
                        else if (regions.size() > 1) {
                            wt_exists = true;
                            wt_regions.add(regions);
                            wt_transcripts_list.add(wt);
                            wt_temp.add(protIdParser(wt.protein_ids));
                        }
                    }
                }

                // If a WT exists for this SV intron → record the event
                if (wt_exists) {
                    int[] exonmm = min_max_exons(wt_regions);
                    int[] basemm = min_max_bases(wt_regions, sv_intron, wt_transcripts_list);

                    int nprots = geneToGetSkips.countUniqueProteins();
                    int ntrans = geneToGetSkips.gene_transcripts.size();

                    String key = geneToGetSkips.gene_id + "_" + sv_intron.toString();
                    if (!skipMap.containsKey(key)) {
                        ExonSkipNotice event = new ExonSkipNotice(
                            geneToGetSkips.gene_id,
                            geneToGetSkips.gene_symbol,
                            geneToGetSkips.chr,
                            geneToGetSkips.strand,
                            nprots,
                            ntrans,
                            sv_intron,
                            wt_regions,
                            wt_temp,    // WT first
                            sv_temp,    // then SV
                            exonmm[0],
                            exonmm[1],
                            basemm[0],
                            basemm[1]
                        );
                        skipMap.put(key, event);
                    }
                }
            }
        }

        exonSkipRegions.addAll(skipMap.values());
        return exonSkipRegions;
    }
 */