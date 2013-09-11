package net.maizegenetics.gwas.modelfitter;

import java.awt.Frame;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

public class StepwiseOLSModelFitterPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(StepwiseOLSModelFitterPlugin.class);

    public StepwiseOLSModelFitterPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

	@Override
	public DataSet performFunction(DataSet input) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ImageIcon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getButtonName() {
		return "Model";
	}

	@Override
	public String getToolTipText() {
		return "Fit multiple markers in a single model (experimental).";
	}

}
