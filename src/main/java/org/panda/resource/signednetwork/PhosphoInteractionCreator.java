package org.panda.resource.signednetwork;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.miner.*;
import org.biopax.paxtools.pattern.util.DifferentialModificationUtil;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class PhosphoInteractionCreator
{
	SIFMiner miner;

	public PhosphoInteractionCreator(SIFMiner miner)
	{
		this.miner = miner;
	}

	public Set<SIFInteraction> create(Match m, IDFetcher fetcher, boolean gainedSite)
	{
		BioPAXElement sourceER = m.get(miner.getSourceLabel(), miner.getPattern());
		BioPAXElement targetER = m.get(miner.getTargetLabel(), miner.getPattern());

		Set<String> sources = fetcher.fetchID(sourceER);
		Set<String> targets = fetcher.fetchID(targetER);

		SIFType sifType = miner.getSIFType();

		Set<SIFInteraction> set = new HashSet<>();

		for (String source : sources)
		{
			for (String target : targets)
			{
				if (source.equals(target)) continue;

				set.add(new SignedSIFInteraction(source, target, sourceER, targetER, sifType,
					new HashSet<>(m.get(((MinerAdapter) miner).getMediatorLabels(), miner.getPattern())),
					new HashSet<>(m.get(((MinerAdapter) miner).getSourcePELabels(), miner.getPattern())),
					new HashSet<>(m.get(((MinerAdapter) miner).getTargetPELabels(), miner.getPattern())),
					DifferentialModificationUtil.collectChangedPhosphorylationSites(
						(PhysicalEntity) m.get("input simple PE", miner.getPattern()),
						(PhysicalEntity) m.get("output simple PE", miner.getPattern()),
						gainedSite)));
			}
		}
		return set;
	}
}
