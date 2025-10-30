public class CDS {
    public final int start;
    public final int end;
    public final String exon_id;

    public CDS(int start, int end, String exon_id){
        this.start = start;
        this.end = end;
        this.exon_id = exon_id;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }
}
