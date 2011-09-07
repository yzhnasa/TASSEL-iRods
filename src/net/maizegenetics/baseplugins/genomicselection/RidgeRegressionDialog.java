package net.maizegenetics.baseplugins.genomicselection;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

public class RidgeRegressionDialog extends JDialog implements ActionListener {
	protected JComboBox cbo1;
	protected JComboBox cbo2;
	protected JTextField txtFrom;
	protected JTextField txtTo;
	protected JTextField txtBy;
	protected boolean runClicked = false;
	protected String firstDataset;
	protected String secondDataset;
	
	private final String CMD_cbo1 = "cbo1";
	private final String CMD_cbo2 = "cbo2";
	private final String CMD_run= "run";
	private final String CMD_quit = "quit";
	public static final String TYPE_PHENO = "phenotypes";
	public static final String TYPE_GENO = "genotypes";

	public static void main(String[] args) {
		RidgeRegressionDialog rrd = new RidgeRegressionDialog("Test Dialog", null, "phenotype data", "genotype data");
		rrd.setVisible(true);
		System.exit(0);
	}
	
	public RidgeRegressionDialog(String title, Frame parentFrame, String firstdataset, String seconddataset) {
		super(parentFrame, title, true);
		firstDataset = firstdataset;
		secondDataset = seconddataset;
		
		setSize(new Dimension(400, 400));
		JPanel panel = new JPanel();
		JScrollPane jsp = new JScrollPane(panel);
		add(jsp);
		panel.setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		
		String[] types = new String[]{TYPE_PHENO, TYPE_GENO};
		cbo1 = new JComboBox(types);
		cbo1.addActionListener(this);
		cbo1.setActionCommand(CMD_cbo1);
		cbo2 = new JComboBox(types);
		cbo2.setSelectedIndex(1);
		cbo2.addActionListener(this);
		cbo2.setActionCommand(CMD_cbo2);
		
		txtFrom = new JTextField("1", 6);
		txtTo = new JTextField("25", 6);
		txtBy = new JTextField("2", 6);

		int yval = 0;
		
		gbc.gridx = 0;
		gbc.gridy = yval++;
		gbc.gridwidth = 1;
		gbc.weightx = 0;
		gbc.weighty = 0;
		gbc.insets = new Insets(10,10,10,10);
		
		panel.add(new JLabel("Dataset Name"), gbc);
		
		gbc.gridx = 1;
		panel.add(new JLabel("Type"), gbc);
		
		gbc.gridy = yval++;
		gbc.gridx = 0;
		gbc.insets = new Insets(10,10,5,10);
		panel.add(new JLabel(firstDataset), gbc);
		gbc.gridx = 1;
		panel.add(cbo1, gbc);
		
		gbc.gridy = yval++;
		gbc.gridx = 0;
		gbc.insets = new Insets(5,10,10,10);
		panel.add(new JLabel(secondDataset), gbc);
		gbc.gridx = 1;
		panel.add(cbo2, gbc);
		
		gbc.gridy = yval++;
		gbc.gridx = 0;
		gbc.gridwidth = 2;
		gbc.insets = new Insets(10,10,10,10);
		//panel.add(new JLabel("Input Parameters"), gbc);

		JPanel parameterPanel = new JPanel(new GridBagLayout());
		parameterPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createLineBorder(Color.BLACK), "Input Parameters"));
		panel.add(parameterPanel, gbc);
		
		gbc.gridy = 0;
		gbc.gridx = 0;
		gbc.gridwidth = 1;
		parameterPanel.add(new JLabel("Vary lambda from"), gbc);
		gbc.gridx = 1;
		parameterPanel.add(txtFrom, gbc);
		
		gbc.gridy = 1;
		gbc.gridx = 0;
		parameterPanel.add(new JLabel("to"), gbc);
		gbc.gridx = 1;
		parameterPanel.add(txtTo, gbc);
		
		gbc.gridy = 2;
		gbc.gridx = 0;
		parameterPanel.add(new JLabel("by"), gbc);
		gbc.gridx = 1;
		parameterPanel.add(txtBy, gbc);
		
		JButton btnRun = new JButton("Run");
		btnRun.addActionListener(this);
		btnRun.setActionCommand(CMD_run);
		JButton btnQuit = new JButton("Quit");
		btnQuit.addActionListener(this);
		btnQuit.setActionCommand(CMD_quit);
		
		gbc.gridy = yval++;
		gbc.gridx = 0;
		panel.add(btnRun, gbc);
		gbc.gridx = 1;
		panel.add(btnQuit, gbc);
		

	}

	public void actionPerformed(ActionEvent event) {
		if (event.getActionCommand().equals(CMD_cbo1)) {
			if (cbo1.getSelectedIndex() == 0) cbo2.setSelectedIndex(1);
			else cbo2.setSelectedIndex(0);
		}
		else if(event.getActionCommand().equals(CMD_cbo2)) {
			if (cbo2.getSelectedIndex() == 0) cbo1.setSelectedIndex(1);
			else cbo1.setSelectedIndex(0);
		}
		else if(event.getActionCommand().equals(CMD_run)) {
			runClicked = true;
			setVisible(false);
		}
		else if(event.getActionCommand().equals(CMD_quit)) {
			setVisible(false);
		}
	}

	public String getTxtFrom() {
		return txtFrom.getText();
	}

	public void setTxtFrom(String from) {
		txtFrom.setText(from);
	}

	public String getTxtTo() {
		return txtTo.getText();
	}

	public void setTxtTo(String to) {
		txtTo.setText(to);
	}

	public String getTxtBy() {
		return txtBy.getText();
	}

	public void setTxtBy(String by) {
		txtBy.setText(by);
	}

	public boolean isRunClicked() {
		return runClicked;
	}

	public String getFirstDataset() {
		return firstDataset;
	}

	public void setFirstDataset(String firstDataset) {
		this.firstDataset = firstDataset;
		
	}
	
}
