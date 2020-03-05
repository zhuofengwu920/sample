

import java.util.ArrayList;
import java.util.Random;




class Mutator extends ArrayList<Changer> {

    // Mutators are lists of Changers

    public void print_mutations(Fragment f) throws Exception {
	int i ;
        Changer ch;
	boolean ready=true;

	StringBuffer buff = new StringBuffer(f.frag);
	for(i = 0; i < size(); i++) {
	    ch =  get(i);
	    ch.make_change(buff);
        }
        System.out.println(">"+f.id+"\t"+f.left+"\t"+f.right+"\n"+buff);
    }


    public void accumulate(Fragment f) throws Exception {
	int i ;
        Changer ch;
	boolean ready=true;;

	for(i = 0; i < size(); i++) {
	    ch =  get(i);
	    ch.accumulate(f);
        }
    }


    public void  global_changes(Mutator localmuts, int num_seqs) 
      throws Exception {
	int i ;
        Changer ch;
	boolean ready=true;;

	for(i = 0; i < size(); i++) {
	    ch =  get(i);
	    ch.global_changes(localmuts, num_seqs);
        }
    }



}





abstract class Changer {

    static Random rnd = new Random();

    protected char dna_letters [] = {'a','c','g','t'};

    public void make_change(StringBuffer buff) throws Exception {
	throw new Exception("Changer.make_change not implemented");
    }

    public void  accumulate(Fragment f) throws Exception {
	throw new Exception("Changer.accumulate not implemented");
    }

    public void  global_changes(Mutator localmuts, int n) throws Exception  {
	throw new Exception("Changer.global_changes not implemented");
    }


    static double get_double(String msg, String val, double a, double b)
	throws Exception {
	double v;
	v = Double.parseDouble(val);
	if (v<a || v>b) {
	    throw new Exception(msg+" out of range. Value: "+v+". Should be between "+
				a+" and "+b);
	}
	return v;
    }


    static int get_int(String msg, String val, int a, int b)
	throws Exception {
	int v;
	v = Integer.parseInt(val);
	if (v<a || v>b) {
	    throw new Exception(msg+" out of range. Should be between "+
				a+" and "+b);
	}
	return v;
    }

 
}


class  BaseChanger extends Changer {


    private int subs, dels, ins, ns;
    private double alpha, beta, gamma;
    private double zeta, xi;

    public BaseChanger(String args[], int curr) throws Exception {
	alpha = get_double("alpha",args[curr],0.0,1.0);
	beta  = get_double("beta",args[curr+1],0.0,80);
	gamma = get_double("gamma",args[curr+2],0.0,1);
	zeta  = get_double("zeta",args[curr+3],0.0,1000);
        xi    = get_double("xi",args[curr+4],0.0,1000);
        subs  = Integer.parseInt(args[curr+5]); // kappa
        dels  = subs+Integer.parseInt(args[curr+6]); // lambda
	ins   = dels+Integer.parseInt(args[curr+7]); // mu
        ns    = ins +Integer.parseInt(args[curr+8]); // nu
    }

    protected StringBuffer base_mute(StringBuffer b, int k) {
	char ch;
        int  r;

	r = rnd.nextInt(ns);
	if (r<subs) {
	    do {
		ch = dna_letters[rnd.nextInt(4)];
            } while (ch==b.charAt(k));
	    b.setCharAt(k,ch);
	} else {
	    if (r<dels) {
	      b.deleteCharAt(k);
	    } else {
		if (r < ins) {
	         ch = dna_letters[rnd.nextInt(4)];
	         b.insert(k,ch);}
	       else 
		   b.setCharAt(k,'N');
	    }
	}
	return b;
    }

    public void make_change(StringBuffer buff) {
	int i, t;
        double cutpt1, cutpt2, cutpt3, cutpt4;
        for(i=0; i<buff.length(); i++) {
            t = buff.length();
	    cutpt1 = alpha;
	    cutpt2 = alpha*(xi-1)*i/t;
            cutpt3 = Math.exp(-Math.sqrt(t-i)/zeta);
            cutpt4 = 0.5*gamma*(1-extmath.tanh(i-beta));
	    
	    if (rnd.nextDouble() < cutpt1 ||
                rnd.nextDouble() < cutpt2 ||
                rnd.nextDouble() < cutpt3 ||
                rnd.nextDouble() < cutpt4)
	       buff = base_mute(buff,i);
	}
    }
}



class Stutter extends Changer {

    double prob;  // Probability that a stutter happens

    // goes through a fragment and adds stutters

    public Stutter (String args [], int curr) throws Exception {
	prob = get_double("eta",args[curr],0.0,100.0);
    }


    private double sqr(double x) {
	return x*x;
    }


    public void make_change(StringBuffer buff) {
	int i, rep=0, stutlen, start;
	char curr, prev='X';
	double cutpt, dice;
	String the_stutter;

	for(i=0; i<buff.length(); i++) {
	    curr = Character.toUpperCase(buff.charAt(i));
	    if (curr != prev) {
		cutpt = 
		    (prev == 'G') ? sqr(0.5 * (1 - Math.cos(rep/prob))): 
		    (prev == 'T') ? sqr(0.5 * (1 - Math.cos(2*rep/prob))):  
		    0;
		dice = rnd.nextDouble();
		if (dice < cutpt) {
		    stutlen = rnd.nextInt(rep*2);
		    start = Math.max(0,i-stutlen);
		    stutlen=i-start+1;
		    the_stutter = buff.substring(start,i);
		    buff.insert(i,the_stutter);
		    i=i+stutlen;
		    rep=0;
		}

		rep = 1;
	    } else 
		rep++;
	    prev=curr;
	}
    }

}


class Ligate extends Changer {

    int n=0;
    ArrayList<Fragment> frags;
    double theta;

    public Ligate(String args[], int index) throws Exception {
	theta = get_double("theta",args[index],0.0,1.0);
	frags = new ArrayList<Fragment>();
    }


    public void accumulate(Fragment f) {
	if (rnd.nextDouble() < theta)
	    frags.add(f);
    }


    public void global_changes(Mutator localmuts, int numseqs) 
      throws Exception {
	int i, j, len;
	Fragment x, y;
	while (frags.size() >= 2) {
	    i=rnd.nextInt(frags.size());
	    x =  frags.get(i);
	    frags.remove(i);
	    j=rnd.nextInt(frags.size());
	    y =  frags.get(j);
	    frags.remove(j);
            x.frag = x.frag + y.frag;
	    len = Math.min(x.frag.length(), 300+rnd.nextInt(200));
            x.frag = x.frag.substring(0,len);
            x.id   = (1+Math.min(x.id,y.id))*numseqs+Math.max(x.id,y.id);
	    localmuts.print_mutations(x);
	} 
    }

}
	    









