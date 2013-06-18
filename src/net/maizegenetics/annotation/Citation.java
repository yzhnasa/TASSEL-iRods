/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.annotation;

/**
 * Provides the publication citation information for an algorithm.  Ultimately, this
 * citation should be passed on to the resulting datasets, so that inventors of cool
 * algoritms get credit for their work.
 * @author edbuckler
 */
public @interface Citation {
    String journalCitation();
    boolean pubEmbargoed();
    String pubEmbargoDate();
}
