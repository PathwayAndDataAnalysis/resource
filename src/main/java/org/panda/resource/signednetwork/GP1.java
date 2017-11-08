package org.panda.resource.signednetwork;

import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.pattern.miner.ControlsStateChangeOfMiner;
import org.biopax.paxtools.pattern.miner.IDFetcher;
import org.biopax.paxtools.pattern.miner.SIFInteraction;

import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class GP1 extends ControlsStateChangeOfMiner
{
	public GP1()
	{
		setType(SignedType.ACTIVATES_GTPASE);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NOT(ConBox.linkToSpecific()), "input PE", "output simple PE");
		p.add(new NOT(ConBox.linkToSpecific()), "output PE", "input simple PE");
		p.add(new Type(Complex.class), "input PE");
		p.add(new Type(Complex.class), "output PE");
		p.add(new OR(
			new MappedConst(
				new AND(
					new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 0),
					new MappedConst(new GTPChangeConstraint(GTPChangeConstraint.Type.GDP_TO_GTP), 1, 2)
				), 0, 1, 2),
			new MappedConst(
				new AND(
					new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 0),
					new MappedConst(new GTPChangeConstraint(GTPChangeConstraint.Type.GTP_TO_GDP), 1, 2)
				), 0, 1, 2)
			), "Control", "input PE", "output PE");
		return p;
	}
}
