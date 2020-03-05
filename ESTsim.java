


import java.io.*;

import java.util.ArrayList;
import java.util.Iterator;




public class ESTsim {

/* ESTsim 
   Scott Hazelhurst
   v0.3    26 August 2008
   v0.2.0  28 January 2003  // Update version number
    0.2.1  29 January 2003 
    0.2.3  19 February 2003
    0.2.4  24 February 2003
    0.2.5  4 March 2003
    0.2.6  27 June 2003
*/

    static String version = "0.3  27 August 2008";


    static Mutator [] create_mutator(String args []) throws Exception {
	String kind;
	int    rate, su, in, de, curr;
	double alpha, beta;
	Mutator mutate [] = new Mutator [2];

        mutate[0] = new Mutator();
        mutate[1] = new Mutator();  
        
        curr = 3;
        while (curr < args.length) {
	    kind = args[curr];
	    if (kind.equals("sbe")) {
		if (curr+9 >= args.length) 
		    throw new Exception("Not enough args for sbe option");
		mutate[0].add(new BaseChanger(args, curr+1));
		curr = curr+10;
		continue;
	    }
	    if (kind.equals("stutter")) {
		if (curr+1 >= args.length)
		    throw new Exception("Not enough args for stutter option");
		mutate[0].add(new Stutter(args, curr+1));
                curr=curr+2;
		continue;
	    }
	    if (kind.equals("ligate")) {
		if (curr+1 >= args.length)
		    throw new Exception("Not enough args for ligate option");
		mutate[1].add(new Ligate(args, curr+1));
		curr=curr+2;
		continue;
	    }
	    throw new Exception("Illegal error model "+kind);
	}
	return mutate;
    }




    public static void main(String args []) throws Exception {

        int       seqno=0, k, num_clones;
	ArrayList<Breaker>  cutters;
        ArrayList<Fragment>  frags;
	SimpleInp baseseqs, cutpoints;
	Mutator   mutate [] = new Mutator [2];
	String    line, fmut;
        Changer   ch;

        if (args[0].equals("-v") || args[0].equals("--version"))
	{ System.out.println(version); 
	    return;}
        

	cutters    = new ArrayList<Breaker>();
	frags      = new ArrayList<Fragment>();
	baseseqs   = new SimpleInp(args[0]);
	cutpoints  = new SimpleInp(args[1]);
	num_clones = Integer.parseInt(args[2]);
        mutate     = create_mutator(args);

    

	do { // read in the cut file
	    line = cutpoints.readLine();
	    if (line == null) break;
	    cutters.add(new Breaker(line));
	} while (true);


        // Produce the local mutations 
	do { // read in the sequences
	    line = baseseqs.readLine();
	    if (line == null) break;
	    if (line.charAt(0) == '>') continue;
	    for(Breaker currcut : cutters ) {
		frags = currcut.breakUp(line,seqno);
		for(Fragment frag : frags ) {
		    mutate[1].accumulate(frag);
		    for(k=0; k<num_clones; k++) 
			mutate[0].print_mutations(frag);
		}
	    }
	    seqno++;
	} while (true);
        // Produce the global mutations
	mutate[1].global_changes(mutate[0], seqno);
    }


}
