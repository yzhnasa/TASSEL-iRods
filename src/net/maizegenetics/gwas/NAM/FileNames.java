package net.maizegenetics.gwas.NAM;

import java.io.File;

public class FileNames {
	public File snps = null;
	public File chrmodel = null;
	public File chrsteps = null;
	public File model = null;
	public File steps = null;
	public File residuals = null;
	public File agpmap = null;
	public File namMarkersByChr = null;
	public File namMarkers = null;
	public File phenotypes = null;
	public File projectedFile = null;
	
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
	public int chromosome = 0;
	
}
