
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.StringTokenizer;



abstract class Cutter  {

    protected int num;
    double rc;
    Random rand = new Random();

    class Counter implements Iterator<Fragment> {
	String seq;
	int base = 1, seqno ;
	public Counter(String seq, int seqno) {
	    this.seq=seq;
	    this.seqno = seqno;
	}
	public boolean hasNext() {
	    return end_seq(seq, base);
	}
	public Fragment next() {
	    int i;
	    Fragment frag;
	    frag = getNext(seq,seqno, base);
	    base = frag.right+1;
	    if (rand.nextDouble() < rc) 
		frag.rc();
	    return frag;
	}

	public void remove(){};
    }


    public  Cutter()  {};

    public  Cutter(String line) {};
    
    abstract Fragment getNext(String seq, int seqno, int from);


    Iterator<Fragment> init(String seq, int seqno) {
        num=0;
	return new Counter(seq, seqno);
    }



    Iterator<Fragment> iterator(String seq, int seqno) {
        num=0;
	return new Counter(seq, seqno);
    }

    protected boolean end_seq(String seq, int base) {
       return base>=0 && base<seq.length();
    }

}


class RFLP extends Cutter {

    String target;

    public RFLP(String init) throws Exception {
	String strinp;
	StringTokenizer tok = new StringTokenizer(init);

	strinp = tok.nextToken();
	if (!strinp.equals("pattern")) 
	    throw new Exception("Illegal data "+strinp);
        strinp = tok.nextToken();
	target = strinp;
	strinp = tok.nextToken();
	rc = Double.parseDouble(strinp);
    }

    Fragment getNext(String seq, int seqno, int from) {
        Fragment f;
	String s;
	int k;
	k = seq.indexOf(target,from);
        s = k<0 ? seq.substring(from) : seq.substring(from,k);
	f = new Fragment(seqno,s,from,from+s.length()-1);
	return f;
    }
}


class RandomCutter extends Cutter  {

    Random rand = new Random();
    String target;
    int    min, max;

    public RandomCutter(String init) throws Exception {
	String strinp;

	StringTokenizer tok = new StringTokenizer(init);

	strinp = tok.nextToken();
	if (!strinp.equals("random")) 
	    throw new Exception("Illegal data "+strinp);
	strinp = tok.nextToken();
	min = Integer.parseInt(strinp);
	strinp = tok.nextToken();
	max = Integer.parseInt(strinp);
	strinp = tok.nextToken();
	rc = Double.parseDouble(strinp);
    }


    Fragment getNext(String seq, int seqno, int from) {
        Fragment f;
        int k;
	int maxlen = 
            Math.min(seq.length()-from-min,max-min);
	k =  (maxlen>0?rand.nextInt(maxlen):maxlen)+min+from;
	if (k+min/3 > seq.length()) k = seq.length();
	f = new Fragment(seqno,seq.substring(from,k),from,k-1);
	return f;
    }

}

class RandomSampleCutter extends Cutter  {

    Random rand = new Random();
    String target;
    int    min, max, cover;

    public RandomSampleCutter(String init) throws Exception {
	String strinp;

	StringTokenizer tok = new StringTokenizer(init);

	strinp = tok.nextToken();
	if (!strinp.equals("samplerandom")) 
	    throw new Exception("Illegal data "+strinp);
	strinp = tok.nextToken();
	min = Integer.parseInt(strinp);
	strinp = tok.nextToken();
	max = Integer.parseInt(strinp);
	strinp = tok.nextToken();
	cover = Integer.parseInt(strinp);
	strinp = tok.nextToken();
	rc = Double.parseDouble(strinp);

    }



     protected boolean end_seq(String seq, int base) {
         return num < 2*cover*seq.length()/(max+min);
    }


    Fragment getNext(String seq, int seqno, int from) {
        Fragment f;
        int k;
	from = Math.max(0,rand.nextInt(seq.length()-min));
	int maxlen = 
            Math.min(seq.length()-from-min,max-min);
	k =  (maxlen>0?rand.nextInt(maxlen):maxlen)+min+from;
	if (k+min/3 > seq.length()) k = seq.length();
	f = new Fragment(seqno,seq.substring(from,k),from,k-1);
        num++;
	return f;
    }

}



class SpliceCutter extends Cutter  {

    ArrayList<Integer> splice_points;
    Random rand = new Random();
    int    fringe; // what area around splice point should not be cut
    int    backward, forward;
    int    min, max;  // range of EST size
    int    cover; // what coverage of cDNA
    int    distrib[];
    int    invert[];


    private boolean near_splice(int x) {
	for (Integer pt : splice_points) 
	    if ((pt-backward <= x) && (x<=pt+forward) ) return true;
	return false;
    }

    public SpliceCutter(String init) throws Exception {
	String strinp, endpoints[];
	int i, j, k, sum=0, last=0;

	splice_points = new ArrayList<Integer>();
	invert = new int [100];
	StringTokenizer tok = new StringTokenizer(init);

        // nsplice sites -1 fringe min max cover rc
	strinp = tok.nextToken();
	if (!strinp.equals("nsplice")) 
	    throw new Exception("Illegal data "+strinp);
	distrib = new int[10];
	do {
	    strinp = tok.nextToken();
	    if (strinp.equals("@")) break;
	    endpoints = strinp.split("\\.\\.");
	    i = Integer.parseInt(endpoints[0]);
	    j = Integer.parseInt(endpoints[1]);
	    last = last+j-i+1;
	    splice_points.add(last);
	} while(true);
	fringe = Integer.parseInt(tok.nextToken());
	backward = fringe/2;
	forward  = fringe-backward;
	min    = Integer.parseInt(tok.nextToken());
	max    = Integer.parseInt(tok.nextToken());
	cover  = Integer.parseInt(tok.nextToken());
	rc     = Double.parseDouble(tok.nextToken());

        for(i=0; i<10; i++) {
	    distrib[i]=Integer.parseInt(tok.nextToken());
	    sum = sum+distrib[i];
	}
	if (sum != 100)
	    throw new Exception("Distribution must equal 100, but is:"+sum);
        k=0;
	sum=0;
	for(i=0; i<10; i++) {
	    sum = sum+distrib[i];
	    while(k<sum){
		invert[k]=i;
		k++;

	    }
	}
    }

     protected boolean end_seq(String seq, int base) {
         return num < 2*cover*seq.length()/(max+min);
    }


    Fragment getNext(String seq, int seqno, int from) {
        Fragment f;
	String frag_s;
        int seqlen,k,  end, r;
	float winlen;

	seqlen = seq.length();
	winlen = (float) seqlen/10+1;
        do {
	    r = rand.nextInt(100);
	    from = invert[r]*((int)winlen)+ rand.nextInt((int)winlen+1);
	    end  = Math.min(seqlen,from+min+rand.nextInt(max-min));
	} while ( (end-from < min*0.5) || near_splice(from) || near_splice(end));

	frag_s = seq.substring(from,end);
	f = new Fragment(seqno,frag_s,from,end-1);
        num++;
	return f;
    }
    
}



class Breaker {
    
    protected Cutter cutter;

    public Breaker(String line) throws Exception {
	switch (line.charAt(0)) {
	case 'r' :  cutter = new RandomCutter(line); break;
        case 's' :  cutter = new RandomSampleCutter(line); break;
	case 'n' :  cutter = new SpliceCutter(line); break;
        default :
          cutter = new RFLP(line);
	}
    }


    ArrayList<Fragment> breakUp(String seq, int seqno) {
	int index=0;

	Fragment curr;
	ArrayList<Fragment> frags = new ArrayList<Fragment>();
        for(Iterator i=cutter.init(seq,seqno); i.hasNext(); ) {
	    curr = (Fragment) i.next();
	    frags.add(curr);
	};
	return frags;
    }


}
