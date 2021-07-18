package org.panda.resource.signednetwork;

import org.biopax.paxtools.pattern.miner.SIFMiner;
import org.biopax.paxtools.pattern.miner.SIFType;

import java.util.Arrays;
import java.util.List;

/**
 * Enumeration of signed and directed SIF relations.
 *
 * @author Ozgun Babur
 */
public enum SignedType implements SIFType
{
	PHOSPHORYLATES("First protein positively affects phosphorylation of the second protein.",
		true, true, true),
	DEPHOSPHORYLATES("First protein negatively affects phosphorylation of the second protein.",
		true, true, false),
	ACETYLATES("First protein positively affects acetylation of the second protein.",
		true, true, true),
	DEACETYLATES("First protein negatively affects acetylation of the second protein.",
		true, true, false),
	METHYLATES("First protein positively affects methylation of the second protein.",
		true, true, true),
	DEMETHYLATES("First protein negatively affects methylation of the second protein.",
		true, true, false),
	UPREGULATES_EXPRESSION("First protein positively affects expression of the second protein.",
		true, false, true),
	DOWNREGULATES_EXPRESSION("First protein negatively affects expression of the second protein.",
		true, false, false),
	ACTIVATES_GTPASE("First protein activates the target GTPase signaling function by one of the several possible " +
		"ways.",
		true, false, true),
	INHIBITS_GTPASE("First protein inhibits the target GTPase signaling function by one of the several possible " +
		"ways.",
		true, false, false),
	PRODUCES("Source protein produces the target metabolite.",
		true, false, true),
	CONSUMES("Source protein consumes the target metabolite.",
		true, false, false),
	USED_TO_PRODUCE("Source metabolite is converted into the target metabolite.",
		true, false, true),
	;

	/**
	 * Constructor with parameters.
	 * @param description description of the edge type
	 * @param directed whether the edge type is directed
	 */
	SignedType(String description, boolean directed, boolean siteSpecific, boolean positive,
		Class<? extends SIFMiner>... miners)
	{
		this.description = description;
		this.directed = directed;
		this.siteSpecific = siteSpecific;
		this.positive = positive;
		this.miners = Arrays.asList(miners);
	}

	/**
	 * Description of the SIF type.
	 */
	private String description;

	/**
	 * Some SIF edges are directed and others are not.
	 */
	private boolean directed;

	/**
	 * Whether this relation is a site-specific event.
	 */
	private boolean siteSpecific;

	/**
	 * We always have a positive and a negative relation of similar types. Like acetylation (positive) and deacetylation
	 * (negative).
	 */
	private boolean positive;

	/**
	 * SIF Miners to use during a search.
	 */
	private List<Class<? extends SIFMiner>> miners;

	/**
	 * Tag of a SIF type is derived from the enum name.
	 * @return tag
	 */
	public String getTag()
	{
		return name().toLowerCase().replaceAll("_", "-");
	}

	/**
	 * Asks if the edge is directed.
	 * @return true if directed
	 */
	public boolean isDirected()
	{
		return directed;
	}

	public boolean isSiteSpecific()
	{
		return siteSpecific;
	}

	public boolean isPositive()
	{
		return this.positive;
	}

	/**
	 * Gets the description of the SIF type.
	 * @return description
	 */
	public String getDescription()
	{
		return description;
	}

	@Override
	public List<Class<? extends SIFMiner>> getMiners()
	{
		return miners;
	}

	public static SignedType typeOf(String tag)
	{
		tag = tag.toUpperCase().replaceAll("-", "_");
		SignedType type = null;
		try
		{
			type = valueOf(tag);
		}
		catch (IllegalArgumentException e){}
		return type;
	}
}
