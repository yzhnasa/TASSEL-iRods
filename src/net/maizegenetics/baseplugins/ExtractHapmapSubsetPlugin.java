package net.maizegenetics.baseplugins;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.TreeSet;
import java.util.regex.Pattern;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;
import org.apache.log4j.xml.DOMConfigurator;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.Utils;


public class ExtractHapmapSubsetPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(ExtractHapmapSubsetPlugin.class);
    private static final Pattern tab = Pattern.compile("\t");

	private int minAlleleCount = 0;
	private String inputFilename = "";
	private String outputFilename = "";
	private String pedFilename = "";
	
	public ExtractHapmapSubsetPlugin(Frame parentFrame) {
		super(parentFrame, false);
	}
	
	@Override
	public DataSet performFunction(DataSet input) {
		extractSubset();
		return null;
	}

	@Override
	public void setParameters(String[] args) {
		if (args == null || args.length < 3) {
			myLogger.error(getUsage());
			return;
		}
		
		int narg = args.length;
		for (int i = 0; i < narg - 1; i++) {
			if (args[i].equals("-h") || args[i].equalsIgnoreCase("-inputfile")) {
				inputFilename = args[++i];
			}
			else if (args[i].equals("-o") || args[i].equalsIgnoreCase("-outputfile")) {
				outputFilename = args[++i];
			}
			else if (args[i].equals("-p") || args[i].equalsIgnoreCase("-pedfile")) {
				pedFilename = args[++i];
			}
			else if (args[i].equals("-a") || args[i].equalsIgnoreCase("-minalleles")) {
				minAlleleCount = Integer.parseInt(args[++i]);
			}
		}
	}

	public static String getUsage() {
		StringBuilder usage = new StringBuilder("The ExtractHapmapSubsetPlugin requires the following parameter:\n");
		usage.append("-h or -inputfile : a hapmap file containing the full data set.\n");
		usage.append("-o or -outputfile : the name of the file that will contain the subset specified by the pedfile.\n");
		usage.append("-p or -pedfile : the pedigree file containing a list of taxa to be extracted.\n");
		usage.append("The following parameter is optional:\n");
		usage.append("-a or -minalleles : only snps with the minimum number of alleles will be extracted. Default = 0.\n\n");
		
		return usage.toString();
	}
	
	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return null;
	}

	@Override
	public String getToolTipText() {
		return null;
	}

	public void extractSubset() {
    	if (inputFilename.length() == 0 || outputFilename.length() == 0 || pedFilename.length() == 0) {
    		myLogger.info("Extract taxa requires three parameters: the input hapmap file, the output hapmap file, and the pedigree file.");
    		myLogger.info(getUsage());
    		System.exit(-1);
    	}
		try {
			LinkedList<String> pedlist = new LinkedList<String>();
			
			System.out.println("Reading taxa file, " + pedFilename);
			BufferedReader br = Utils.getBufferedReader(pedFilename);
			String input;
			while ((input = br.readLine()) != null) {
				pedlist.add(input);
			}
			
			br.close();
			
			br = Utils.getBufferedReader(inputFilename);
			System.out.println("Reading hapmap file, " + inputFilename);
			input = br.readLine();
			String[] info = tab.split(input);
			int n = info.length;
			LinkedList<Integer> colList = new LinkedList<Integer>();
			Collections.sort(pedlist);
			
			for (int i = 11; i < n; i++) {
				if (Collections.binarySearch(pedlist, info[i]) > -1) colList.add(i);
			}
			
			BufferedWriter bw = Utils.getBufferedWriter(outputFilename);
			bw.write(info[0]);
			for (int i = 1; i < 11; i++) {
				bw.write("\t");
				bw.write(info[i]);
			}
			for (Integer col:colList) {
				bw.write("\t");
				bw.write(info[col]);
			}
			bw.newLine();
			
			if (minAlleleCount > 0) {
				while ((input = br.readLine()) != null) {
					info = tab.split(input);
					TreeSet<String> alleles = new TreeSet<String>();
					StringBuilder sb = new StringBuilder(info[0]);
					for (int i = 1; i < 11; i++) sb.append("\t").append(info[i]);
					for (Integer col:colList) {
						alleles.add(info[col]);
						sb.append("\t").append(info[col]);
					}
					int nAlleles = alleles.size();
					boolean hasN = alleles.contains("N");
					if ((hasN && nAlleles >= minAlleleCount + 1) || (!hasN && nAlleles >= minAlleleCount)) {
						bw.write(sb.toString());
						bw.newLine();
					}
				}
			} else {
				while ((input = br.readLine()) != null) {
					info = tab.split(input);
					StringBuilder sb = new StringBuilder(info[0]);
					for (int i = 1; i < 11; i++) sb.append("\t").append(info[i]);
					for (Integer col:colList) {
						sb.append("\t").append(info[col]);
					}
					bw.write(sb.toString());
					bw.newLine();
				}
			}
			
			br.close();
			bw.close();
	
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		System.out.println("Finished. Output written to " + outputFilename);
	}
	
}
