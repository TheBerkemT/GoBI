public class IntronRegion {
    int start;
    int end;
    
    public IntronRegion(int s, int e){
        this.start = s;
        this.end = e;
    }

    public int getEnd() {
        return end;
    }
    public int getStart() {
        return start;
    }

    @Override
    public String toString() {
        return ""+this.start+":"+this.end;
    }
}
