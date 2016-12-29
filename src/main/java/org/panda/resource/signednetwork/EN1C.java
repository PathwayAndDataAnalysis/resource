package org.panda.resource.signednetwork;

import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.NOT;
import org.biopax.paxtools.pattern.constraint.Type;
import org.biopax.paxtools.pattern.miner.ControlsExpressionMiner;

/**
 * Negative expression pattern 1 with a complex.
 *
 * @author Ozgun Babur
 */
public class EN1C extends EN1
{
	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NOT(new Type(Complex.class)), "TF PE");
		return p;
	}
}
