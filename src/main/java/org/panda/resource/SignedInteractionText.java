package org.panda.resource;

import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.CollectionUtil;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class SignedInteractionText
{
	String source;
	String target;
	SignedType type;
	String mediators;
	Set<String> sites;
	Set<String> origSites;
	Set<String> remSites;

	public SignedInteractionText(String line)
	{
		String[] t = line.split("\t");
		source = t[0];
		type = SignedType.typeOf(t[1]);
		target = t[2];
		if (t.length > 3 && !t[3].isEmpty())
		{
			mediators = t[3];
		}
		if (t.length > 4 && !t[4].isEmpty())
		{
			sites = new HashSet<>(Arrays.asList(t[4].split(";")));
			origSites = new HashSet<>(sites);
		}
	}

	public String key()
	{
		return source + " " + type.getTag() + " " + target;
	}

	public String signlessKey()
	{
		return source +
			(type == SignedType.PHOSPHORYLATES || type == SignedType.DEPHOSPHORYLATES ?
				" phospho " : " expression ") +
			target;
	}

	public void subtractSites(SignedInteractionText inter)
	{
		if (inter.sites != null && sites != null)
		{
			remSites = CollectionUtil.getIntersection(sites, inter.sites);
			sites.removeAll(remSites);
		}
	}

	public void addSites(SignedInteractionText sit)
	{
		if (sites == null)
		{
			if (sit.sites != null) sites = new HashSet<>(sit.sites);
		}
		else if (sit.sites != null) sites.addAll(sit.sites);
	}

	public void addMediators(SignedInteractionText inter)
	{
		if (mediators == null)
		{
			mediators = inter.mediators;
		}
		else if (inter.mediators != null)
		{
			mediators = CollectionUtil.merge(CollectionUtil.getUnion(
				Arrays.asList(mediators.split(";")), Arrays.asList(inter.mediators.split(";"))), ";");
		}
	}

	@Override
	public String toString()
	{
		String s = source + "\t" + type.getTag() + "\t" + target;
		if (mediators != null) s += "\t" + mediators;
		if (sites != null)
		{
			if (mediators == null)
			{
				s += "\t";
			}

			s += "\t" + CollectionUtil.merge(sites, ";");
		}
		return s;
	}

	public String toStringWOType()
	{
		String s = source + "\t" + target;
		if (mediators != null) s += "\t" + mediators;
		if (sites != null)
		{
			if (mediators == null)
			{
				s += "\t";
			}

			s += "\t" + CollectionUtil.merge(sites, ";");
		}
		return s;
	}

	public Set<String> getSites()
	{
		return sites;
	}

	public void setSites(Set<String> sites)
	{
		this.sites = sites;
	}

	public Set<String> getOrigSites()
	{
		return origSites;
	}

	public void setOrigSites(Set<String> origSites)
	{
		this.origSites = origSites;
	}

	public Set<String> getRemSites()
	{
		return remSites;
	}

	public void setRemSites(Set<String> remSites)
	{
		this.remSites = remSites;
	}

	public SignedType getType()
	{
		return type;
	}

	public void merge(SignedInteractionText sit)
	{
		assert source.equals(sit.source);
		assert target.equals(sit.target);
		assert type.equals(sit.type);

		addMediators(sit);
		addSites(sit);
	}
}
