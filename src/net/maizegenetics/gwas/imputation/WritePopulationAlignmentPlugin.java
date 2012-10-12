package net.maizegenetics.gwas.imputation;

import java.awt.Frame;
import java.io.File;
import java.util.Arrays;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

public class WritePopulationAlignmentPlugin extends AbstractPlugin {
	private static final Logger myLogger = Logger.getLogger(WritePopulationAlignmentPlugin.class);
	boolean mergeAlignments = true;
	boolean writeParentCalls = true;
	boolean outputDiploid = false;
	String baseFileName;
	
	public WritePopulationAlignmentPlugin(Frame parentFrame) {
		super(parentFrame, false);
	}
	
	@Override
	public DataSet performFunction(DataSet input) {
		List<Datum> theData = input.getDataOfType(PopulationData.class);
		if (mergeAlignments) {
			String filename = baseFileName + ".hmp.txt";
			Alignment[] allOfTheAlignments = new Alignment[theData.size()];
			int count = 0;
			for (Datum datum:theData) { 
				PopulationData family = (PopulationData) datum.getData();
				int nsnps = (int) family.snpIndex.cardinality(); //only output the snps set in family.snpIndex
				int[] snpsToOutput = new int[nsnps];
				int snpcount = 0;
				for (int s = 0; s < nsnps; s++) {
					if (family.snpIndex.fastGet(s)) snpsToOutput[snpcount++] = s;
				}
				snpsToOutput = Arrays.copyOf(snpsToOutput, snpcount);
				Alignment outputAlignment;
				if (writeParentCalls) {
					outputAlignment = FilterAlignment.getInstance(family.imputed, snpsToOutput);
				}
				else {
					outputAlignment = FilterAlignment.getInstance(family.original, snpsToOutput);
				}

				allOfTheAlignments[count++] = outputAlignment;
			}
			Alignment alignment = MutableSingleEncodeAlignment.getInstance(allOfTheAlignments);
			ExportUtils.writeToHapmap(alignment, outputDiploid, filename, '\t', null);
		} else {
			for (Datum datum:theData) {
				PopulationData family = (PopulationData) datum.getData();
				int nsnps = (int) family.snpIndex.cardinality();
				int[] snpsToOutput = new int[nsnps];
				int snpcount = 0;
				for (int s = 0; s < nsnps; s++) {
					if (family.snpIndex.fastGet(s)) snpsToOutput[snpcount++] = s;
				}
				snpsToOutput = Arrays.copyOf(snpsToOutput, snpcount);
				Alignment outputAlignment;
				String filename = baseFileName + ".family." + family.name + ".hmp.txt";
				if (writeParentCalls) {
					outputAlignment = FilterAlignment.getInstance(family.imputed, snpsToOutput);
				}
				else {
					outputAlignment = FilterAlignment.getInstance(family.original, snpsToOutput);
				}
				ExportUtils.writeToHapmap(outputAlignment, outputDiploid, filename, '\t', null);
			}
		}
		
		fireDataSetReturned(new PluginEvent(input, input.getCreator().getClass()));
		return input;
	}

	@Override
	public void setParameters(String[] args) {
		if (args == null || args.length == 0) {
			myLogger.error(getUsage());
			return;
		}
		
		int narg = args.length;
		for (int i = 0; i < narg - 1; i++) {
			if (args[i].equals("-f") || args[i].equalsIgnoreCase("-file")) {
				baseFileName = args[++i];
			}
			else if (args[i].equals("-m") || args[i].equalsIgnoreCase("-merge")) {
				String val = args[++i];
				if (val.toUpperCase().startsWith("T")) mergeAlignments = true;
				else mergeAlignments = false;
			}
			else if (args[i].equals("-p") || args[i].equalsIgnoreCase("-parentCalls")) {
				String val = args[++i];
				if (val.toUpperCase().startsWith("T")) writeParentCalls = true;
				else writeParentCalls = false;
			}
			else if (args[i].equals("-d") || args[i].equalsIgnoreCase("-diploid")) {
				String val = args[++i];
				if (val.toUpperCase().startsWith("T")) outputDiploid = true;
				else writeParentCalls = false;
			}
			else if (args[i].equals("?")) myLogger.error(getUsage());
		}
	}

	public void setMergeAlignments(boolean mergeAlignments) {
		this.mergeAlignments = mergeAlignments;
	}

	public void setWriteParentCalls(boolean writeParentCalls) {
		this.writeParentCalls = writeParentCalls;
	}

	public void setOutputDiploid(boolean outputDiploid) {
		this.outputDiploid = outputDiploid;
	}

	public void setBaseFileName(String baseFileName) {
		this.baseFileName = baseFileName;
	}

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return "Write Populations";
	}

	@Override
	public String getToolTipText() {
		return null;
	}

	private String getUsage() {
		StringBuilder usage = new StringBuilder("The WritePopulationAlignmentPlugin requires the following parameter:\n");
		usage.append("-f or -file : The base file name for the ouput. .hmp.txt will be appended.\n");
		usage.append("The following parameters are optional:\n");
		usage.append("-m or -merge : if true families are merged into a single file, if false each family is output to a separate file (default = true)");
		usage.append("-p or -parentCalls : if true output is A/C/M where parent 1 is A and parent 2 is M, if false outputs nucleotides (default = true)");
		usage.append("-d or -diploid : if true output is AA/CC/AC, if false output is A/C/M");
		usage.append("? : print the parameter list.");

		return usage.toString();
	}

}
