package net.maizegenetics.gwas.NAMgwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Pattern;

public class FileNames {
	public File[] snps = new File[10];
	public File[] chrmodel = new File[10];
	public File[] chrsteps = new File[10];
	public File model = null;
	public File steps = null;
	public File residuals = null;
	public File agpmap = null;
	public File[] namMarkersByChr = new File[10];
	public File namMarkers = null;
	public File phenotypes = null;
	public File[] projectedFile = new File[10];
	
	public double enterlimit = Double.NaN;
	public double exitlimit = Double.NaN;
	public double[] limits = null;
	public int traitnumber = -1;
	public int iterations = 0;
	public String replacement;
	public String analysis = "";
	public int maxsnps = 100;
	public boolean threaded = true;
	public boolean permute = false;
	public boolean bootstrapPermutation = false;
	public boolean randomizeSnpOrder = false;
	public int startIteration = 0;
	
	public Pattern ws;
	public Pattern equal;
	public Pattern quote;
	
	public static final String TYPE_SNP = "snp";
	public static final String TYPE_RESIDUAL = "residual";
	public static final String TYPE_AGPMAP = "agpmap";
	public static final String TYPE_MODEL = "model";
	public static final String TYPE_STEPS = "steps";
	public static final String TYPE_NAMGENO = "namgeno";
	public static final String TYPE_PHENOTYPE = "phenotype";
	public static final String TYPE_PROJECTED = "projected";
	public static final String TYPE_RANDOM = "random";
	
	public static final String NAME_ENTERLIMIT = "enterlimit";
	public static final String NAME_LIMITBYCHR = "enterlimits";
	public static final String NAME_TRAITNUMBER = "traitnumber";
	public static final String NAME_SNPLIST = "snplist";
	public static final String NAME_ITERATIONS = "iterations";
	public static final String NAME_WITH_REPLACEMENT = "replacement";
	public static final String NAME_ANALYSIS = "analysis";
	public static final String NAME_MAXSNPS = "maxsnps";
	public static final String NAME_THREADED = "threaded";
	public static final String NAME_PERMUTE = "permute";
	public static final String NAME_START_ITERATION = "startiteration";
	public static final String NAME_SNP_ORDER = "snporder";
	
	public static final String ATTR_TYPE = "type";
	public static final String ATTR_NAME = "name";
	public static final String ATTR_CHR = "chr";
	public static final String ATTR_VALUE = "value";
	
	public FileNames(String inputFile) {
		//create patterns
		ws = Pattern.compile("\\s");
		equal = Pattern.compile("=");
		quote = Pattern.compile("\"");

		BufferedReader br;
		
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String input= br.readLine();
			while (input != null) {
				if (!parseFile(input)) parseParameter(input);
				input = br.readLine();
			}
			br.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	public FileNames() {}
	
	public boolean parseFile(String inline) {
		String type = null;
		String name = null;
		int chr = -1;
		
		String[] splitstr = ws.split(inline, 2);
		if (!splitstr[0].equalsIgnoreCase("<File")) return false;
		splitstr[1] = splitstr[1].trim();
		
		if (splitstr[1].endsWith(">")) {
			int n = splitstr[1].length();
			splitstr[1] = splitstr[1].substring(0, n - 1);
		}
		
		String attr;
		String value;
		while (splitstr.length > 1 && splitstr[1].length() > 0) {
			splitstr = equal.split(splitstr[1], 2);
			attr = splitstr[0].trim();
			splitstr[1] = splitstr[1].trim();
			if (splitstr[1].startsWith("\"")) {
				splitstr = quote.split(splitstr[1].substring(1), 2);
				splitstr[1] = splitstr[1].trim();
			}
			else splitstr = ws.split(splitstr[1], 2);
			value = splitstr[0];
			
			if (attr.equalsIgnoreCase(ATTR_CHR)) {chr = Integer.parseInt(value);}
			else if (attr.equalsIgnoreCase(ATTR_NAME)) {name = value;}
			else if (attr.equalsIgnoreCase(ATTR_TYPE)) {type = value;}
		}
		
		if (type.equalsIgnoreCase(TYPE_AGPMAP)) {agpmap = new File(name);}
		else if (type.equalsIgnoreCase(TYPE_MODEL)) {
			if (chr == -1) model = new File(name);
			else chrmodel[chr - 1] = new File(name);
		}
		else if (type.equalsIgnoreCase(TYPE_NAMGENO)) {
			if (chr == -1) namMarkers = new File(name);
			namMarkersByChr[chr -1 ] = new File(name);
		}
		else if (type.equalsIgnoreCase(TYPE_RESIDUAL)) {residuals = new File(name);}
		else if (type.equalsIgnoreCase(TYPE_SNP)) {snps[chr - 1] = new File(name);}
		else if (type.equalsIgnoreCase(TYPE_STEPS)) {
			if (chr == -1) steps = new File(name);
			else chrsteps[chr - 1] = new File(name);
		}
		else if (type.equalsIgnoreCase(TYPE_PHENOTYPE)) {phenotypes = new File(name);}
		else if (type.equalsIgnoreCase(TYPE_PROJECTED)) {projectedFile[chr - 1] = new File(name); }
		else System.err.println("In file list: unknown type for name=" + name);

		return true;
	}
	
	public boolean parseParameter(String inline) {
		String name = "";
		String paramvalue = "";
		String type = "";
		
		String[] splitstr = ws.split(inline, 2);
		if (!splitstr[0].equalsIgnoreCase("<Parameter")) return false;
		splitstr[1] = splitstr[1].trim();
		
		if (splitstr[1].endsWith(">")) {
			int n = splitstr[1].length();
			splitstr[1] = splitstr[1].substring(0, n - 1);
		}
		
		String attr;
		String value;
		while (splitstr.length > 1 && splitstr[1].length() > 0) {
			splitstr = equal.split(splitstr[1], 2);
			attr = splitstr[0].trim();
			splitstr[1] = splitstr[1].trim();
			splitstr = ws.split(splitstr[1], 2);
			value = splitstr[0];
			
			if (attr.equalsIgnoreCase(ATTR_NAME)) {name = value;}
			else if (attr.equalsIgnoreCase(ATTR_VALUE)) {paramvalue = value;}
			else if (attr.equalsIgnoreCase(ATTR_TYPE)) {type = value;}
		}

		if (name.equalsIgnoreCase(NAME_ENTERLIMIT)) enterlimit = Double.parseDouble(paramvalue);
		else if (name.equalsIgnoreCase(NAME_TRAITNUMBER)) traitnumber = Integer.parseInt(paramvalue);
		else if (name.equalsIgnoreCase(NAME_ITERATIONS)) iterations = Integer.parseInt(paramvalue);
		else if (name.equalsIgnoreCase(NAME_WITH_REPLACEMENT)) replacement = paramvalue;
		else if (name.equalsIgnoreCase(NAME_ANALYSIS)) {
			analysis = paramvalue;
			if (type.toLowerCase().startsWith("boot")) bootstrapPermutation = true;
		}
		else if (name.equalsIgnoreCase(NAME_MAXSNPS)) maxsnps = Integer.parseInt(paramvalue);
		else if (name.equalsIgnoreCase(NAME_THREADED) && paramvalue.toLowerCase().startsWith("f")) threaded = false;
		else if (name.equalsIgnoreCase(NAME_SNP_ORDER) && paramvalue.toLowerCase().startsWith("r")) randomizeSnpOrder = true;
		else if (name.equalsIgnoreCase(NAME_PERMUTE) && paramvalue.toLowerCase().startsWith("t")) {
			permute = true;
		}
		else if (name.equalsIgnoreCase(NAME_LIMITBYCHR)) {
			String[] limitstr = paramvalue.split(",");
			int n = limitstr.length;
			limits = new double[n];
			for (int i = 0; i < n; i++) {
				try {
					limits[i] = Double.parseDouble(limitstr[i]);
				} catch (Exception e) {
					limits[i] = 1e-6;
				}
				System.out.println("Enter limit " + i + " = " + limits[i]);
			}
		}
		else if (name.equalsIgnoreCase(NAME_START_ITERATION)) startIteration = Integer.parseInt(paramvalue);
		
		return true;
	}
	
}
