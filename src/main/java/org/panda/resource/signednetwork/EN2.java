package org.panda.resource.signednetwork;

import org.biopax.paxtools.pattern.Pattern;

/**
 * @author Ozgun Babur
 */
public class EN2 extends EP2
{
	public EN2()
	{
		setType(SignedType.DOWNREGULATES_EXPRESSION);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.removeLastConstraint();
		p.add(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), "Control");
		return p;
	}
}
