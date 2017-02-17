package org.panda.resource.signednetwork;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.pattern.miner.*;
import org.biopax.paxtools.pattern.util.DifferentialModificationUtil;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class PP4 extends CSCOThroughControllingSmallMoleculeMiner
{
	public PP4()
	{
		setType(SignedType.PHOSPHORYLATES);
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = super.constructPattern();
		p.add(new NOT(ConBox.linkToSpecific()), "input PE", "output simple PE");
		p.add(new NOT(ConBox.linkToSpecific()), "output PE", "input simple PE");
		p.add(new OR(
			new MappedConst(
				new AND(
					new MappedConst(
						new OR(
							new MappedConst(
								new AND(
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 0),
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 1)
								), 0, 1),
							new MappedConst(
								new AND(
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 0),
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 1)
								), 0, 1)
						), 0, 1
					),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.GAIN, "phospho"), 2, 3)
				), 0, 1, 2, 3),
			new MappedConst(
				new AND(
					new MappedConst(
						new OR(
							new MappedConst(
								new AND(
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 0),
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 1)
								), 0, 1),
							new MappedConst(
								new AND(
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 0),
									new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 1)
								), 0, 1)
						), 0, 1
					),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.LOSS, "phospho"), 2, 3)
				), 0, 1, 2, 3)
			), "upper Control", "Control", "input simple PE", "output simple PE");
		return p;
	}

	@Override
	public Set<SIFInteraction> createSIFInteraction(Match m, IDFetcher fetcher)
	{
		PhosphoInteractionCreator pic = new PhosphoInteractionCreator(this);
		return pic.create(m, fetcher, useGainedSite(m));
	}

	protected boolean useGainedSite(Match m)
	{
		boolean ctrlSign1 = new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE)
			.satisfies(m, getPattern().indexOf("upper Control"));

		boolean ctrlSign2 = new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE)
			.satisfies(m, getPattern().indexOf("Control"));

		boolean ctrlSign = ctrlSign1 == ctrlSign2;

		boolean edgeSign = getSIFType().equals(SignedType.PHOSPHORYLATES);

		return edgeSign == ctrlSign;
	}
}
