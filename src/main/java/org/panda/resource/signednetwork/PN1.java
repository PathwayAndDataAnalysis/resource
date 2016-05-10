package org.panda.resource.signednetwork;

import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.pattern.miner.ControlsStateChangeOfMiner;

/**
 * @author Ozgun Babur
 */
public class PN1 extends PP1
{
	public PN1()
	{
		setType(SignedType.DEPHOSPHORYLATES);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.removeLastConstraint();
		p.add(new OR(
			new MappedConst(
				new AND(
					new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 0),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.LOSS, "phospho"), 1, 2)
				), 0, 1, 2),
			new MappedConst(
				new AND(
					new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 0),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.GAIN, "phospho"), 1, 2)
				), 0, 1, 2)
			), "Control", "input simple PE", "output simple PE");
		return p;
	}
}
