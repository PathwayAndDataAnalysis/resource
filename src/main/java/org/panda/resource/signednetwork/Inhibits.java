package org.panda.resource.signednetwork;

import org.biopax.paxtools.pattern.miner.SIFMiner;
import org.biopax.paxtools.pattern.miner.SIFType;

import java.util.List;

/**
 * @author Ozgun Babur
 */
public class Inhibits implements SIFType
{
	@Override
	public String getTag()
	{
		return "inhibits";
	}

	@Override
	public boolean isDirected()
	{
		return true;
	}

	@Override
	public String getDescription()
	{
		return "Source protein binds and inhibits the target";
	}

	@Override
	public List<Class<? extends SIFMiner>> getMiners()
	{
		throw new UnsupportedOperationException("This SIF type should not be asked for miners.");
	}

	@Override
	public int hashCode()
	{
		return getTag().hashCode();
	}

	@Override
	public boolean equals(Object obj)
	{
		return obj instanceof Inhibits && ((Inhibits) obj).getTag().equals(getTag());
	}
}
