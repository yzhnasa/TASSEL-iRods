package net.maizegenetics.analysis.distance;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextArea;

public class KinshipMethodDialog extends JDialog {
	JRadioButton radioEndelman;
	JRadioButton radioIBS;
	boolean isCancelled = true;
	
	public KinshipMethodDialog(Frame parent) {
		super(parent, "Choose Kinship Method", true);
		
		radioEndelman = new JRadioButton("Scaled IBS", true);
		radioIBS = new JRadioButton("Pairwise IBS", false);
		ButtonGroup bg = new ButtonGroup();
		bg.add(radioEndelman);
		bg.add(radioIBS);
		
        JPanel panel = new JPanel();
        panel.setLayout(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();

        JButton okButton = new JButton();
        okButton.setActionCommand("Ok");
        okButton.setText("Ok");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	isCancelled = false;
                setVisible(false);
            }
        });
        
        JButton cancelButton = new JButton();
        cancelButton.setText("CancelClose");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	isCancelled = true;
                setVisible(false);
            }
        });
        
        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.X_AXIS));
        buttonPanel.add(okButton);
        buttonPanel.add(Box.createHorizontalStrut(15));
        buttonPanel.add(cancelButton);
        
        StringBuilder helpText = new StringBuilder("The Scaled IBS method produces a kinship matrix that is scaled to give a reasonable estimate of additive genetic variance.");
        helpText.append(" The Pairwise IBS method, which is the method used by TASSEL version 4, may result in an inflated estimate of genetic variance.");
        helpText.append(" Either will do a good job of controlling population structure in MLM.");
        final JTextArea txtInfo = new JTextArea(helpText.toString());
        txtInfo.setLineWrap(true);
        txtInfo.setWrapStyleWord(true);
        txtInfo.setPreferredSize(new Dimension(500, 150));
        txtInfo.setEditable(false);

        gbc.gridx = 0;
        gbc.gridy = 0;
        panel.add(Box.createVerticalStrut(15));
        
        gbc.gridy++;
        gbc.anchor = GridBagConstraints.CENTER;
        gbc.insets = new Insets(5,0,5,0);
        panel.add(radioEndelman, gbc);
        gbc.gridy++;
        panel.add(radioIBS, gbc);

        gbc.gridy++;
        gbc.gridx = 0;
        gbc.anchor = GridBagConstraints.CENTER;
        gbc.insets = new Insets(5,5,5,5);
        panel.add(buttonPanel, gbc);
        
        gbc.gridy++;
        gbc.gridx = 0;
        gbc.insets = new Insets(20,0,10,0);
        panel.add(txtInfo, gbc);
        
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        
        getContentPane().add(panel, BorderLayout.CENTER);
        pack();
        setVisible(true);
	}
	
	public boolean isEndelmanSelected() {
		return radioEndelman.isSelected();
	}
	
	public boolean isCancelled() {
		return isCancelled;
	}
}
