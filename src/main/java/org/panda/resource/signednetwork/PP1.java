package org.panda.resource.signednetwork;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.pattern.MappedConst;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.pattern.miner.*;
import org.biopax.paxtools.pattern.util.DifferentialModificationUtil;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class PP1 extends ControlsStateChangeOfMiner
{
	public PP1()
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
					new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.POSITIVE), 0),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.GAIN, "phospho"), 1, 2)
				), 0, 1, 2),
			new MappedConst(
				new AND(
					new MappedConst(new ControlSignConstraint(ControlSignConstraint.Sign.NEGATIVE), 0),
					new MappedConst(new ModificationChangeConstraint(ModificationChangeConstraint.Type.LOSS, "phospho"), 1, 2)
				), 0, 1, 2)
			), "Control", "input simple PE", "output simple PE");
		return p;
	}

	@Override
	public Set<SIFInteraction> createSIFInteraction(Match m, IDFetcher fetcher)
	{
		BioPAXElement sourceER = m.get(this.getSourceLabel(), getPattern());
		BioPAXElement targetER = m.get(this.getTargetLabel(), getPattern());

		Set<String> sources = fetcher.fetchID(sourceER);
		Set<String> targets = fetcher.fetchID(targetER);

		SIFType sifType = this.getSIFType();

		Set<SIFInteraction> set = new HashSet<>();

		for (String source : sources)
		{
			for (String target : targets)
			{
				if (source.equals(target)) continue;

				set.add(new SignedSIFInteraction(source, target, sourceER, targetER, sifType,
					new HashSet<>(m.get(getMediatorLabels(), getPattern())),
					new HashSet<>(m.get(getSourcePELabels(), getPattern())),
					new HashSet<>(m.get(getTargetPELabels(), getPattern())),
					DifferentialModificationUtil.collectChangedPhosphorylationSites(
						(PhysicalEntity) m.get("input simple PE", getPattern()),
						(PhysicalEntity) m.get("output simple PE", getPattern()),
						getSIFType().equals(SignedType.PHOSPHORYLATES))));
			}
		}
		return set;
	}
}
