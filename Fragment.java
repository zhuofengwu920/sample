



class Fragment {
    int id;
    public String frag;
    int left, right;
    boolean rcd;

    public Fragment(int identity,String fragment,int lft,int rgt) {
	id   = identity;
	frag = fragment;
	left = lft;
	right= rgt;
	rcd = false;
    }


    static char compl(char d) {
	switch(d) {
	case 'a' : return 't';
	case 'c' : return 'g';
	case 't' : return 'a';
	case 'g' : return 't';
	default : return d;
	}
    }


    public void rc() {
       char ch;
       int i,len = frag.length();
       StringBuffer tmp = new StringBuffer(frag);
       for(i=0; i<len/2; i++) {
          ch = tmp.charAt(i);
	  tmp.setCharAt(i,compl(tmp.charAt(len-1-i)));
	  tmp.setCharAt(len-1-i, compl(ch));
       }
       rcd = true;
       frag = new String(tmp);
    }
}





