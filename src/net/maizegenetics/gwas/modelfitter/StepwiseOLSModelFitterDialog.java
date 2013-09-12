package net.maizegenetics.gwas.modelfitter;

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JRootPane;

import net.maizegenetics.pal.alignment.Alignment;

public class StepwiseOLSModelFitterDialog extends JDialog implements ActionListener{
	JCheckBox jboxNested = new JCheckBox("Nest markers within main factor", null, false);
	JList<String> listMainEffects = null;
	JLabel lblchoose = new JLabel("Choose a main effect to nest within:");
	JButton btnOk = new JButton("btnOK");
	JButton btnCancel = new JButton("btnCancel");
	String[] mainEffects = null;
	int indexOfSelectedEffect = 0;
	boolean wasCancelled = false;
	
	public StepwiseOLSModelFitterDialog(String[] mainEffects, Frame parentFrame) {
		super();
		this.mainEffects = mainEffects;
		this.setTitle("Choose Nested Model");
		this.setModalityType(ModalityType.APPLICATION_MODAL);
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);
        Container contentPane = getContentPane();
//        BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
        JPanel myPanel = new JPanel();
        myPanel.setLayout(new BoxLayout(myPanel, BoxLayout.Y_AXIS));
        myPanel.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));
        contentPane.add(myPanel);
        
        
		int numberOfMainEffects = mainEffects.length;
		if (numberOfMainEffects > 1) {
			myPanel.add(Box.createVerticalStrut(20));
			StringBuilder question = new StringBuilder("The dataset contains more than one main effect.");
			JLabel lblquestion = new JLabel(question.toString());
			lblquestion.setAlignmentX(CENTER_ALIGNMENT);
			myPanel.add(lblquestion);
			myPanel.add(Box.createVerticalStrut(20));
			jboxNested.setAlignmentX(CENTER_ALIGNMENT);
			myPanel.add(jboxNested);
			jboxNested.addActionListener(this);
			jboxNested.setActionCommand("nest");
			listMainEffects = new JList<String>(mainEffects);
			lblchoose.setVisible(true);
			listMainEffects.setVisible(true);
			listMainEffects.setBorder(BorderFactory.createLineBorder(Color.BLACK));
			myPanel.add(Box.createVerticalStrut(30));
			myPanel.add(lblchoose);
			lblchoose.setAlignmentX(Component.CENTER_ALIGNMENT);
			myPanel.add(Box.createVerticalStrut(10));
			myPanel.add(listMainEffects);
			listMainEffects.setAlignmentX(Component.CENTER_ALIGNMENT);
			
			myPanel.add(Box.createVerticalStrut(30));
			myPanel.add(btnOk);
			btnOk.setAlignmentX(CENTER_ALIGNMENT);
			btnOk.addActionListener(this);
			btnOk.setActionCommand("ok");
			
			myPanel.add(btnCancel);
			btnCancel.setAlignmentX(CENTER_ALIGNMENT);
			btnCancel.addActionListener(this);
			btnCancel.setActionCommand("cancel");
			
			
		} else {
			StringBuilder question = new StringBuilder("The dataset contains only one main effect, ");
			question.append(mainEffects[0]);
			myPanel.add(new JLabel(question.toString()));
			myPanel.add(Box.createVerticalStrut(20));
			myPanel.add(jboxNested);
		}
		
		pack();
		lblchoose.setVisible(false);
		listMainEffects.setVisible(false);
		setVisible(true);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Action Command: " + e.getActionCommand());
		if (e.getActionCommand().equals("nest")) {
			if (jboxNested.isSelected()) {
				lblchoose.setVisible(true);
				listMainEffects.setVisible(true);
			} else {
				lblchoose.setVisible(false);
				listMainEffects.setVisible(false);
			}
		} else if (e.getActionCommand().equals("ok")) {
			setVisible(false);
		} else if (e.getActionCommand().equals("cancel")) {
			wasCancelled = true;
			setVisible(false);
		}
	
	}
	
	public boolean isNested() {
		return jboxNested.isSelected();
	}
	
	public String getNestedEffect() {
		if (mainEffects == null || mainEffects.length == 0) return null;
		if (mainEffects.length == 1) return mainEffects[0];
		return listMainEffects.getSelectedValue();
	}
	
	public boolean wasCancelled() {
		return wasCancelled;
	}
}
