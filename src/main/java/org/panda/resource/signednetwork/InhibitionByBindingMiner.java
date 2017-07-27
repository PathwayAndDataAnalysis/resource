package org.panda.resource.signednetwork;

import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.pattern.Match;
import org.biopax.paxtools.pattern.Pattern;
import org.biopax.paxtools.pattern.constraint.*;
import org.biopax.paxtools.pattern.miner.AbstractSIFMiner;

import java.util.Set;

import static org.biopax.paxtools.pattern.constraint.ConBox.*;

/**
 * @author Ozgun Babur
 */
public class InhibitionByBindingMiner extends AbstractSIFMiner
{
	ControlToControllerER cc;

	public InhibitionByBindingMiner()
	{
		super(new Inhibits(), "inhibits", "Source protein inhibits activity of the target protein by binding.");
		cc = new ControlToControllerER();
	}

	@Override
	public Pattern constructPattern()
	{
		Pattern p = new Pattern(Interaction.class, "Process");
		p.add(interToControl(), "Process", "Control1");
		p.add(interToControl(), "Process", "Control2");
		p.add(equal(false), "Control1", "Control2");
		p.add(new IsSubsetAndReverse(), "Control1", "Control2");
		p.add(new ControlToControllerER(), "Control2", "ER Source");
		p.add(new NOT(new ControlToControllerER()), "Control1", "ER Source");
		p.add(new ControlToControllerER(), "Control1", "ER Target");
		p.add(new ControlToControllerER(), "Control2", "ER Target");
		p.add(equal(false), "ER Source", "ER Target");
		return p;
	}

	@Override
	public String getSourceLabel()
	{
		return "ER Source";
	}

	@Override
	public String getTargetLabel()
	{
		return "ER Target";
	}

	@Override
	public String[] getMediatorLabels()
	{
		return new String[]{"Control1", "Control2"};
	}

	class IsSubsetAndReverse extends ConstraintAdapter
	{
		public IsSubsetAndReverse()
		{
			super(2);
		}

		@Override
		public boolean satisfies(Match match, int... ind)
		{
			Control c1 = (Control) match.get(ind[0]);
			Control c2 = (Control) match.get(ind[1]);

			if (getSign(c1) * getSign(c2) != -1) return false;

			Set<EntityReference> e1 = cc.getRelatedERs(c1);
			Set<EntityReference> e2 = cc.getRelatedERs(c2);

			return e2.containsAll(e1) && e2.size() > e1.size();
		}

		private int getSign(Control cont)
		{
			ControlType type = cont.getControlType();
			if (type == null) return 0;
			if (type.toString().startsWith("A")) return 1;
			return -1;
		}
	}

}
