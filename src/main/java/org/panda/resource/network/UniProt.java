package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;
import org.panda.utility.Kronometre;
import org.panda.utility.TermCounter;
import org.panda.utility.graph.SiteSpecificGraph;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Serves the PhosphoNetworks database.
 *
 * @author Ozgun Babur
 */
public class UniProt extends FileServer
{
	static UniProt instance;

	static SiteSpecificGraph graph;

	public static UniProt get()
	{
		if (instance == null) instance = new UniProt();
		return instance;
	}

	public SiteSpecificGraph getGraph()
	{
		return graph;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"UniProt-network.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		graph = new SiteSpecificGraph("UniProt", SignedType.PHOSPHORYLATES.getTag());

		Scanner scanner = new Scanner(new File(locateInBase(getLocalFilenames()[0])));
		String target = null;
		while (scanner.hasNextLine())
		{
			String line = scanner.nextLine();

			if (line.startsWith(">"))
			{
				target = line.substring(1);
			}
			else
			{
				String[] token = line.split("\t");
				String source = token[2];
				String site = token[1];

				if (!graph.hasRelation(source, target))
				{
					graph.putRelation(source, target, "", site);
				}
				else
				{
					graph.addSite(source, target, site);
				}
			}
		}

		return true;
	}

	public static void convertFromUniProtXML() throws ParserConfigurationException, SAXException, IOException
	{
		Kronometre k = new Kronometre();
		String inFile = "/media/babur/6TB1/uniprot_sprot.xml";
		String outFile = "/home/babur/Downloads/UniProt-network.txt";

		SAXParserFactory factory = SAXParserFactory.newInstance();
		SAXParser saxParser = factory.newSAXParser();
		UPHandler handler = new UPHandler();
		saxParser.parse(inFile, handler);
		handler.tc.print();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		handler.graph.write(writer);
		writer.close();
		k.stop();
		k.print();
	}

	static class UPHandler extends DefaultHandler
	{

		SiteSpecificGraph graph = new SiteSpecificGraph("UniProt", SignedType.PHOSPHORYLATES.getTag());

		String[] source = null;
		String target = null;
		String modification = null;
		Set<String> refHolders = new HashSet<>();
		Map<String, String> pmIDMap = new HashMap<>();

		List<Relation> relList = new ArrayList<>();

		boolean inGene = false;
		boolean willReadGeneName = false;
		boolean inFeature = false;
		boolean inLocation = false;

		TermCounter tc = new TermCounter();

		@Override
		public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException
		{
			if (qName.equals("entry"))
			{
				target = null;
				pmIDMap.clear();
				relList.clear();
			}
			else if (qName.equals("gene"))
			{
				inGene = true;
			}
			else if (inGene &&
				qName.equals("name") && attributes.getValue("type") != null &&
				attributes.getValue("type").equals("primary"))
			{
				willReadGeneName = true;
			}
			else if (qName.equals("feature") && attributes.getValue("type") != null)
			{
				String type = attributes.getValue("type");
				if (type.equals("modified residue"))
				{
					String description = attributes.getValue("description");
					if (description != null && description.contains("; by "))
					{
						modification = description.substring(0, description.indexOf(";"));
						tc.addTerm(modification);
						source = description.substring(description.indexOf("; by ") + 5).split(", ");

						String evidence = attributes.getValue("evidence");
						if (evidence != null && !evidence.isEmpty()) Collections.addAll(refHolders, evidence.split(" "));

						inFeature = true;
					}
				}
			}
			else if (qName.equals("location") && inFeature)
			{
				inLocation = true;
			}
			else if (inFeature && inLocation && qName.equals("position") && attributes.getValue("position") != null)
			{
				String res = attributes.getValue("position");
				if (modification.equals("Phosphoserine")) res = "S" + res;
				else if (modification.equals("Phosphothreonine")) res = "T" + res;
				else if (modification.equals("Phosphotyrosine")) res = "Y" + res;
				else if (modification.equals("Phosphohistidine")) res = "H" + res;

				for (String s : source)
				{
					Relation rel = new Relation();
					rel.source = s;
					rel.target = target;
					rel.residue = res;
					rel.pubmed.addAll(refHolders);
					relList.add(rel);
				}

			}
		}

		@Override
		public void endElement(String uri, String localName, String qName) throws SAXException
		{
			if (qName.equals("entry"))
			{
				for (Relation rel : relList)
				{
					for (String s : new ArrayList<>(rel.pubmed))
					{
						rel.pubmed.remove(s);
						if (pmIDMap.containsKey(s)) rel.pubmed.add("Pubmed:" + s);
					}

					graph.putRelation(rel.source, rel.target, CollectionUtil.merge(rel.pubmed, ";"), rel.residue);
				}
			}
			else if (qName.equals("gene"))
			{
				inGene = false;
			}
			else if (qName.equals("feature"))
			{
				inFeature = false;
				refHolders.clear();
			}
			else if (qName.equals("location"))
			{
				inLocation = false;
			}
			willReadGeneName = false;
		}

		@Override
		public void characters(char[] ch, int start, int length) throws SAXException
		{
			if (willReadGeneName)
			{
				target = getString(ch, start, length);
			}
		}

		private String getString(char[] ch, int start, int length)
		{
			String s = new String(ch);
			s = s.substring(start, start + length);
			return s;
		}
	}

	static class Relation
	{
		String source;
		String target;
		String residue;
		Set<String> pubmed = new HashSet<>();
	}

	public static void main(String[] args) throws IOException, SAXException, ParserConfigurationException
	{
		convertFromUniProtXML();
//		SiteSpecificGraph pn = get().getGraph();
//		boolean contains = pn.getDownstream("MAPK14").contains("MAPKAPK2");
//		System.out.println("contains = " + contains);
	}
}
