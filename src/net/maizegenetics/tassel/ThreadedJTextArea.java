/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
package net.maizegenetics.tassel;

import javax.swing.*;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */

public class ThreadedJTextArea extends JTextArea implements Runnable {
  private String newText;
  private Thread theThread;

  public ThreadedJTextArea() {
   // theThread = new Thread(this);
  }

  public void start(String nt) {
    newText=nt;
    theThread = new Thread(this);
    theThread.start();
  }
  public void run() {
    setText(newText);
//    repaint();
    theThread=null;
    /**@todo: Implement this java.lang.Runnable method*/
  //  throw new java.lang.UnsupportedOperationException("Method run() not yet implemented.");
  }
}
