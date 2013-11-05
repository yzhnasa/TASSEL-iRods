package net.maizegenetics.tassel;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;

import net.maizegenetics.prefs.TasselPrefs;

/**
 * @author terryc
 */
public class PreferencesDialog extends JDialog {
    
    private final static int TEXT_FIELD_WIDTH = 10;
    private final static Font HEADING_FONT = new Font(null, Font.BOLD, 14);
    private ButtonGroup myButtonGroup = new ButtonGroup();
    private JTextField myMaxRetainAlleles = new JTextField(TEXT_FIELD_WIDTH);
    private JCheckBox myRetainRareAlleles = new JCheckBox("Retain Rare Alleles");
    
    public PreferencesDialog() {
        super((Frame) null, "Preferences...", true);
        try {
            init();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private void init() throws Exception {
        
        setTitle("Preferences...");
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);
        
        Container contentPane = getContentPane();
        
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.LEFT_ALIGNMENT);
        
        result.add(Box.createRigidArea(new Dimension(1, 30)));
        
        result.add(getGenotypeStoringPanel());
        
        result.add(Box.createRigidArea(new Dimension(1, 10)));
        
        result.add(getButtons());
        
        result.add(Box.createRigidArea(new Dimension(1, 20)));
        
        contentPane.add(result);
        
        pack();
        
        setResizable(false);
        
    }
    
    private JPanel getGenotypeStoringPanel() {
        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        
        JLabel alignPrefsLabel = new JLabel("Alignment Preferences...");
        alignPrefsLabel.setFont(HEADING_FONT);
        result.add(alignPrefsLabel);
        result.add(Box.createRigidArea(new Dimension(1, 10)));
        
        myMaxRetainAlleles.setText(String.valueOf(TasselPrefs.getAlignmentMaxAllelesToRetain()));
        result.add(getLine("Max Alleles Retained: ", myMaxRetainAlleles));
        
        myRetainRareAlleles.setSelected(TasselPrefs.getAlignmentRetainRareAlleles());
        result.add(myRetainRareAlleles);
        
        return result;
    }
    
    private JPanel getLine(String label, JTextField ref) {
        
        JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));
        
        result.add(new JLabel(label));
        ref.setEditable(true);
        ref.setHorizontalAlignment(JTextField.LEFT);
        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
        ref.setMaximumSize(ref.getPreferredSize());
        result.add(ref);
        
        return result;
        
    }
    
    private JPanel getButtons() {
        
        JButton okButton = new JButton();
        JButton cancelButton = new JButton();
        
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });
        
        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });
        
        JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));
        
        result.add(okButton);
        
        result.add(cancelButton);
        
        return result;
        
    }
    
    private void okButton_actionPerformed(ActionEvent e) {
        
        try {
            int maxNumAlleles = Integer.valueOf(myMaxRetainAlleles.getText());
            if ((maxNumAlleles < 1) || (maxNumAlleles > 14)) {
                JOptionPane.showMessageDialog(this, "Max Alleles Retained must be between 1 and 14 inclusive", "Error", JOptionPane.ERROR_MESSAGE);
                return;
            }
            TasselPrefs.putAlignmentMaxAllelesToRetain(maxNumAlleles);
        } catch (NumberFormatException ex1) {
            JOptionPane.showMessageDialog(this, "Max Alleles Retained must be integer number", "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        
        TasselPrefs.putAlignmentRetainRareAlleles(myRetainRareAlleles.isSelected());
        
        setVisible(false);
    }
    
    private void cancelButton_actionPerformed(ActionEvent e) {
        
        myMaxRetainAlleles.setText(String.valueOf(TasselPrefs.getAlignmentMaxAllelesToRetain()));
        myRetainRareAlleles.setSelected(TasselPrefs.getAlignmentRetainRareAlleles());
        
        setVisible(false);
    }
}
