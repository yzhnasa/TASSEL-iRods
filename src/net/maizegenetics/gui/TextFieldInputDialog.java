package net.maizegenetics.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;


public class TextFieldInputDialog extends JDialog {
	JLabel[] labels;
	JTextField[] textfields;
	JButton btnRun;
	JButton btnCancel;
	String message;
	boolean runClicked = false;
	
	public static void main(String[] args) {
		GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
		String[] fonts = ge.getAvailableFontFamilyNames();
		JLabel[] lab = new JLabel[fonts.length];
		JTextField[] text = new JTextField[fonts.length];
		for (int i = 0; i < fonts.length; i++) {
			lab[i] = new JLabel(fonts[i]);
			lab[i].setFont(new Font(fonts[i], Font.PLAIN, 14));
			text[i] = new JTextField(10); 
		}
		System.out.println(lab[0].getFont().getFontName() + ", size = " + lab[0].getFont().getSize());
		TextFieldInputDialog cd = new TextFieldInputDialog(lab,text, "the message");
		cd.setTitle("Font family names");
		cd.setVisible(true);
		System.out.println("text box 1: " + text[0].getText());
		System.out.println("text box 2: " + text[1].getText());
		System.out.println("runClicked = " + cd.wasRunClicked());
		System.exit(0);
	}
	
	public TextFieldInputDialog(JLabel[] labels, JTextField[] textfields, String message) {
		super();
		this.labels = labels;
		this.textfields = textfields;
		this.message = message;
		init();
	}
	
	public TextFieldInputDialog() {
		labels = new JLabel[3];
		textfields = new JTextField[3];
		
		for (int i = 0; i < 3; i++) {
			labels[i] = new JLabel();
			textfields[i] = new JTextField(10);
		}
		
		labels[0].setText("Label 1");
		labels[1].setText("Label 2");
		labels[2].setText("Label 3");
		message = "The first column of the data set must be the analysis variable. " +
        "The second column must be the population number. " +
        "The remainder of the columns must be markers. " +
        "Do you want to include population in the model as a separate term?";
		message = null;
		init();
	}
	
	public void init() {
		setModal(true);
		setTitle("Custom Dialog");
		int numberOfItems = labels.length;
		setSize(new Dimension(400,400));
		JPanel panel = new JPanel();
		JScrollPane jsp = new JScrollPane(panel);
		panel.setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.gridheight = 1;
		gbc.gridwidth = 1;
		gbc.anchor = GridBagConstraints.SOUTH;
		gbc.fill = GridBagConstraints.BOTH;
		gbc.weightx = 0;
		gbc.weighty = 0;
		gbc.insets = new Insets(10,10,10,10);
		gbc.ipadx = 0;
		gbc.ipady = 0;

		JButton btnRun = new JButton("Run");
		JButton btnCancel = new JButton("Cancel");
		
		btnRun.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent arg0) {
//				JOptionPane.showMessageDialog(null, "run button pressed");
				runClicked = true;
				setVisible(false);
			}
			
		});
		
		btnCancel.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				runClicked = false;
				setVisible(false);
			}
		});

		int firsty = 0;
		
		if (message != null && message.length() > 0) {
		    gbc.gridwidth=2;
		    JTextArea jat = new JTextArea(message);
		    jat.setLineWrap(true);
		    jat.setWrapStyleWord(true);
		    jat.setBackground(labels[0].getBackground());
		    panel.add(jat, gbc);
		    firsty++;
		    gbc.gridwidth=1;
		}
		for (int i = 0; i < numberOfItems; i++) {
			gbc.gridx = 0;
			gbc.gridy = i + firsty;
			panel.add(labels[i], gbc);
			gbc.gridx = 1;
			panel.add(textfields[i], gbc);
		}
		
		gbc.gridx = 0;
		gbc.gridy = numberOfItems + firsty;
		gbc.insets = new Insets(20, 5, 20, 15);
		panel.add(btnRun, gbc);
		
		gbc.gridx = 1;
		gbc.insets = new Insets(20, 15, 20, 5);
		panel.add(btnCancel, gbc);
		
		panel.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));
		add(jsp);
		
	}

	public boolean wasRunClicked() {
		return runClicked;
	}
	
	
}
