package org.panda.resource.signednetwork;

import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.constraint.ConstraintAdapter;
import org.biopax.paxtools.pattern.util.DifferentialModificationUtil;

import java.util.HashSet;
import java.util.Set;

/**
 * This class checks if the two given complexes has a GTP GDP change event.
 *
 * var0: Substrate Complex
 * Var1: Product Complex
 *
 * @author Ozgun Babur
 */
public class GTPChangeConstraint extends ConstraintAdapter
{
	/**
	 * Gain or loss?
	 */
	protected Type type;

	/**
	 * Constructor with the desired change.
	 *
	 * @param type either GTP to GDP, or GDP to GTP
	 */
	public GTPChangeConstraint(Type type)
	{
		super(2);

		this.type = type;
	}

	/**
	 * Checks the complexes if they have differential GTP and GDP members, in the desired configuration.
	 * @param match current pattern match
	 * @param ind mapped indices
	 * @return true if GTP and GDP exists as desired
	 */
	@Override
	public boolean satisfies(Match match, int... ind)
	{
		Complex c1 = (Complex) match.get(ind[0]);
		Complex c2 = (Complex) match.get(ind[1]);

		Set<String> n1 = collectSmallMolNames(c1);
		Set<String> n2 = collectSmallMolNames(c2);

		if (type == Type.GDP_TO_GTP)
		{
			return n1.contains("GDP") && !n1.contains("GTP") && n2.contains("GTP") && !n2.contains("GDP");
		}
		else if (type == Type.GTP_TO_GDP)
		{
			return n1.contains("GTP") && !n1.contains("GDP") && n2.contains("GDP") && !n2.contains("GTP");
		}
		else
		{
			throw new RuntimeException("Invalid type = " + type);
		}
	}

	/**
	 * Gets the all the names of the small molecules in the given complex.
	 * @param c the complex
	 * @return names of member small molecules
	 */
	private Set<String> collectSmallMolNames(Complex c)
	{
		Set<String> names = new HashSet<>();

		for (SimplePhysicalEntity spe : c.getSimpleMembers())
		{
			if (spe instanceof SmallMolecule)
			{
				names.addAll(spe.getName());

				EntityReference er = spe.getEntityReference();

				if (er instanceof SmallMoleculeReference)
				{
					names.addAll(er.getName());
				}
			}
		}
		return names;
	}

	/**
	 * Possible configuration types.
	 */
	public enum Type
	{
		GDP_TO_GTP,
		GTP_TO_GDP
	}
}
