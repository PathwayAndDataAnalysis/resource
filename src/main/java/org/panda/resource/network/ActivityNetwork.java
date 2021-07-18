package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.resource.siteeffect.Feature;
import org.panda.utility.CollectionUtil;
import org.panda.utility.TermCounter;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.SiteSpecificGraph;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

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
		return new String[]{"SPIKE.sif"};
	}

	@Override
	public boolean load() throws IOException
	{
		posGraph = new DirectedGraph("Activation", "activates");
		negGraph = new DirectedGraph("Inhibition", "inhibits");

		// Load SPIKE

//		posGraph.load(locateInBase(getLocalFilenames()[0]), Collections.singleton("activates"));
//		negGraph.load(locateInBase(getLocalFilenames()[0]), Collections.singleton("inhibits"));

		// Inferring from SignedPC phospho phospGraph

		SiteSpecificGraph posP = new SiteSpecificGraph("Phosphorylation", SignedType.PHOSPHORYLATES.getTag());

		posP.merge(SignedPC.get().getGraph(SignedType.PHOSPHORYLATES));
		posP.merge(IPTMNet.get().getGraph(SignedType.PHOSPHORYLATES));
		posP.merge(PhosphoNetworks.get().getGraph());

		SiteSpecificGraph negP = (SiteSpecificGraph) SignedPC.get().getGraph(SignedType.DEPHOSPHORYLATES);

		SiteEffectCollective sec = new SiteEffectCollective();

		harvest(posP, sec, posGraph, negGraph);
		harvest(negP, sec, negGraph, posGraph);



		return true;
	}


	private void harvest(SiteSpecificGraph phosG, SiteEffectCollective sec, DirectedGraph posEffGraph, DirectedGraph negEffGraph)
	{
		phosG.getOneSideSymbols(true).forEach(s -> phosG.getDownstream(s).forEach(t ->
			phosG.getSites(s, t).forEach(site ->
			{
				int e = getEffect(site, t, sec);
				if (e == 1) posEffGraph.putRelation(s, t, phosG.getMediatorsInString(s, t));
				else if (e == -1) negEffGraph.putRelation(s, t, phosG.getMediatorsInString(s, t));
			})));
	}

	private int getEffect(String site, String gene, SiteEffectCollective sec)
	{
		Integer e = sec.getEffect(gene, site, Feature.PHOSPHORYLATION);
		if (e != null && e != 0) return e;

		String aa = site.substring(0, 1);
		int loc = Integer.valueOf(site.substring(1));

		e = sec.getEffect(gene, aa + (loc + 1), Feature.PHOSPHORYLATION);
		if (e != null && e != 0) return e;

		e = sec.getEffect(gene, aa + (loc - 1), Feature.PHOSPHORYLATION);
		if (e != null && e != 0) return e;

		if (aa.equals("T") || aa.equals("S"))
		{
			aa = aa.equals("T") ? "S" : "T";

			e = sec.getEffect(gene, aa + (loc), Feature.PHOSPHORYLATION);
			if (e != null && e != 0) return e;

			e = sec.getEffect(gene, aa + (loc + 1), Feature.PHOSPHORYLATION);
			if (e != null && e != 0) return e;

			e = sec.getEffect(gene, aa + (loc - 1), Feature.PHOSPHORYLATION);
			if (e != null && e != 0) return e;
		}

		return 0;
	}

	private static void parseSPIKE() throws IOException, SAXException, ParserConfigurationException
	{
		DirectedGraph posGraph = new DirectedGraph("Activation", "activates");
		DirectedGraph negGraph = new DirectedGraph("Inhibition", "inhibits");

		TermCounter tc = new TermCounter();

		Map<String, String> idToGene = new HashMap<>();
		Map<String, Set<String>> groupToMembers = new HashMap<>();


		DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
		Document doc = dBuilder.parse("SPIKE.xml");
		doc.getDocumentElement().normalize();

		NodeList list = doc.getElementsByTagName("Gene");
		for (int i = 0; i < list.getLength(); i++)
		{
			Node item = list.item(i);
			NamedNodeMap attributes = item.getAttributes();
			String name = attributes.getNamedItem("name").getNodeValue();
			String id = attributes.getNamedItem("id").getNodeValue();
			idToGene.put(id, name);
		}

		list = doc.getElementsByTagName("Group");
		for (int i = 0; i < list.getLength(); i++)
		{
			Node item = list.item(i);
			NamedNodeMap attributes = item.getAttributes();
			String id = attributes.getNamedItem("id").getNodeValue();
			groupToMembers.put(id, new HashSet<>());

			NodeList children = item.getChildNodes();
			for (int j = 0; j < children.getLength(); j++)
			{
				Node child = children.item(j);
				if (child.getNodeName().equals("Member"))
				{
					String memID = child.getAttributes().getNamedItem("ref").getNodeValue();
					groupToMembers.get(id).add(memID);
				}
			}
		}

		list = doc.getElementsByTagName("Regulation");
		for (int i = 0; i < list.getLength(); i++)
		{
			Node item = list.item(i);

			String effect = item.getAttributes().getNamedItem("effect").getNodeValue();
			if (effect.equals("3")) continue;
			assert effect.equals("1") || effect.equals("2") : "A new effect?: " + effect;

			tc.addTerm(effect);

			String sourceID = null;
			String targetID = null;
			Set<String> pmIDs = new HashSet<>();

			NodeList children = item.getChildNodes();
			for (int j = 0; j < children.getLength(); j++)
			{
				Node child = children.item(j);
				if (child.getNodeName().equals("Source"))
				{
					sourceID = child.getAttributes().getNamedItem("ref").getNodeValue();
				}
				else if (child.getNodeName().equals("PhysicalTarget"))
				{
					targetID = child.getAttributes().getNamedItem("ref").getNodeValue();
				}
				else if (child.getNodeName().equals("Reference"))
				{
					pmIDs.add("PMID:" + child.getAttributes().getNamedItem("pmid").getNodeValue());
				}
			}

			for (String source : getNames(sourceID, idToGene, groupToMembers))
			{
				for (String target : getNames(targetID, idToGene, groupToMembers))
				{
					DirectedGraph graph = effect.equals("1") ? posGraph : negGraph;
					graph.putRelation(source, target, CollectionUtil.merge(pmIDs, ";"));
				}
			}
		}

		tc.print();

		posGraph.printStats();
		negGraph.printStats();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get("SPIKE.sif"));
		posGraph.write(writer);
		negGraph.write(writer);
		writer.close();
	}

	private static Set<String> getNames(String id, Map<String, String> idToGene, Map<String, Set<String>> groupToMembers)
	{
		if (idToGene.containsKey(id)) return Collections.singleton(idToGene.get(id));

		if (groupToMembers.containsKey(id))
		{
			Set<String> set = new HashSet<>();
			for (String memID : groupToMembers.get(id))
			{
				set.addAll(getNames(memID, idToGene, groupToMembers));
			}
			return set;
		}

		return Collections.emptySet();
	}


	public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException
	{
		DirectedGraph graph = get().getPositiveGraph();
		graph.write("activity-pos.sif");
		graph = get().getNegativeGraph();
		graph.write("activity-neg.sif");
	}
}
