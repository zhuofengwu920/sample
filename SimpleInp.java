


//* SimpleInp puts a wrapper on the 


import java.awt.*;
import java.io.*;
import java.lang.*;

/**
   Provides a simple wrapper to SystemInput and other input classes to allow
   simple input from the console and from files.
   To input from the console use
   -- d = SimpleInp.consoleReadInt();
   -- s = SimpleInp.consoleReadLine();
   -- c = SimpleInp.consoleReadChar();
   to read in an integer, string and char respectively
   
   To input from a file,
   1. need to first declare a reference to SimpleInp object e.g.
       SimpleInp myfile;
   2. create a new object of the type, giving as the argument the name
      of the file on disk
      e.g. myfile = new SimpleInp("rainfall.dat");
   3. Read from the file using the methods readInt, readChar or readLine.
*/
   


    public class SimpleInp   {

	private FileInputStream instr;
	private BufferedReader dinp;
        private InputStreamReader dreader;   




	public static BufferedReader
                 Console =  
                    new BufferedReader(new InputStreamReader(System.in));

	public static String consoleReadLine() throws IOException {
	    return Console.readLine();
	}
   
   
	public static int consoleReadInt() throws IOException {
	    String data;
   	       
	    data = Console.readLine();
   	       
	    return Integer.parseInt(data);
	}
   
	public static char consoleReadChar () throws IOException {
	    int d;
	    d = Console.read();
	    return (char) d;
	}	       

	public  SimpleInp (String name)  throws IOException {
	    instr = new FileInputStream(name);
            dreader = new InputStreamReader(instr);
	    dinp  = new BufferedReader(dreader);
	}
   
	public void finalize () throws IOException {
	    instr.close();
	}
    
	public  String readLine() throws IOException {
	    return dinp.readLine();
	}
   
	public int readInt() throws IOException {
	    String data;
   	       
	    data = dinp.readLine();
   	       
	    return Integer.parseInt(data);
	}
   
   
	public char readChar () throws IOException {
	    int d;
	    d = dinp.read();
	    return (char) d;
	}	       


    }

