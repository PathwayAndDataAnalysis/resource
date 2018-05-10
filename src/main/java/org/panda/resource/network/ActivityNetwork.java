package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import java.util.Set;

/**
 * Serves the PhosphoNetworks database.
 *
 * @author Ozgun Babur
 */
public class ActivityNetwork extends FileServer
{
	static ActivityNetwork instance;

	static DirectedGraph posGraph;
	static DirectedGraph negGraph;

	public static ActivityNetwork get()
	{
		if (instance == null) instance = new ActivityNetwork();
		return instance;
	}

	public DirectedGraph getPositiveGraph()
	{
		return posGraph;
	}

	public DirectedGraph getNegativeGraph()
	{
		return negGraph;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{};
	}

	@Override
	public boolean load() throws IOException
	{
		posGraph = new DirectedGraph("Activation", "activates");
		negGraph = new DirectedGraph("Inhibition", "inhibits");

		PhosphoGraph posP = new PhosphoGraph("Phosphorylation", SignedType.PHOSPHORYLATES.getTag());

		posP.merge(SignedPC.get().getGraph(SignedType.PHOSPHORYLATES));
		posP.merge(IPTMNet.get().getGraph());
		posP.merge(PhosphoNetworks.get().getGraph());

		PhosphoGraph negP = (PhosphoGraph) SignedPC.get().getGraph(SignedType.DEPHOSPHORYLATES);

		SiteEffectCollective sec = new SiteEffectCollective();

		harvest(posP, sec, posGraph, negGraph);
		harvest(negP, sec, negGraph, posGraph);

		return true;
	}

	private void harvest(PhosphoGraph phosG, SiteEffectCollective sec, DirectedGraph posEffGraph, DirectedGraph negEffGraph)
	{
		phosG.getOneSideSymbols(true).forEach(s -> phosG.getDownstream(s).forEach(t ->
			phosG.getSites(s, t).forEach(site ->
			{
				int e = getEffect(site, t, sec);
				if (e == 1) posEffGraph.putRelation(s, t);
				else if (e == -1) negEffGraph.putRelation(s, t);
			})));
	}

	private int getEffect(String site, String gene, SiteEffectCollective sec)
	{
		Integer e = sec.getEffect(gene, site);
		if (e != null && e != 0) return e;

		String aa = site.substring(0, 1);
		int loc = Integer.valueOf(site.substring(1));

		e = sec.getEffect(gene, aa + (loc + 1));
		if (e != null && e != 0) return e;

		e = sec.getEffect(gene, aa + (loc - 1));
		if (e != null && e != 0) return e;

		if (aa.equals("T") || aa.equals("S"))
		{
			aa = aa.equals("T") ? "S" : "T";

			e = sec.getEffect(gene, aa + (loc));
			if (e != null && e != 0) return e;

			e = sec.getEffect(gene, aa + (loc + 1));
			if (e != null && e != 0) return e;

			e = sec.getEffect(gene, aa + (loc - 1));
			if (e != null && e != 0) return e;
		}

		return 0;
	}

	public static void main(String[] args)
	{
		Set<String> inh = ActivityNetwork.get().getNegativeGraph().getDownstream("EGFR");
		System.out.println("inh = " + inh);
	}
}
