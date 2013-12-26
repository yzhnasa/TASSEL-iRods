package net.maizegenetics.phenotype;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.regex.Pattern;

public class ImportPhenotypeUtils {
	private static final Pattern sep = Pattern.compile("\\s+");
	private static final Pattern equal = Pattern.compile("=");
	private static final Pattern brack = Pattern.compile("[<>]");
	
	
	public static Phenotype importPhenotypeFile(String filename) throws IOException {
		String[] traitnames = null;
		LinkedList<String[]> headerList = new LinkedList<String[]>();
		String[] use = null;
		boolean alldata = false;
		boolean allfactor = false;
		boolean allcovar = false;
		
		BufferedReader br = new BufferedReader(new FileReader(filename));
		
		//process header rows and count the non-blank rows
		int numberOfDataLines = 0;
		String inputline = br.readLine().trim();
		while (inputline != null) {
			inputline = inputline.trim();
			String[] parsedline;
			if (inputline.startsWith("<")) {
				int close = inputline.indexOf('>');
				String tag = inputline.substring(1, close);
				if (tag.toUpperCase().startsWith("TRAIT")) {
					parsedline = sep.split(inputline);
					int n = parsedline.length - 1;
					traitnames = new String[n];
					for (int c = 0; c < n; c++) traitnames[c] = parsedline[c + 1];
				} else if (tag.toUpperCase().startsWith("HEADER")) {
					parsedline = sep.split(inputline);
					int n = parsedline.length - 1;
					String[] head = new String[n];
					for (int c = 0; c < n; c++) head[c] = parsedline[c + 1];
					int equalpos = head[0].indexOf('=');
					close = head[0].indexOf('>');
					head[0] = head[0].substring(equalpos + 1, close);
				} else if (tag.toUpperCase().startsWith("USE")) {
					parsedline = sep.split(inputline);
					int n = parsedline.length - 1;
					use = new String[n];
					for (int c = 0; c < n; c++) use[c] = parsedline[c + 1];
				} else if (tag.toUpperCase().startsWith("DATA")) {
					alldata = true;
				} else if (tag.toUpperCase().startsWith("COVAR")) {
					allcovar = true;
				} else if (tag.toUpperCase().startsWith("FACTOR")) {
					allfactor = true;
				}
			} else if (!inputline.startsWith("#")) {
				numberOfDataLines++;
			}
			
			inputline = br.readLine().trim();
		}
		br.close();
		
		//reopen the file and read the data lines
		br = new BufferedReader(new FileReader(filename));
		while ((inputline = br.readLine()) != null) {
			inputline = inputline.trim();
			if (inputline.length() > 0 && !inputline.startsWith("#") && !inputline.startsWith("<")) {
				//process a line of data
				String[] parsedline = sep.split(inputline);
				
			}
		}
		
		br.close();
		return null;
	}
}
