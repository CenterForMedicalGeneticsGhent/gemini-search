import re

REGION_REGEX = re.compile(r"^([^:]+):(\d+)-(\d+)$")


def parse_region(region_str):
    """Parse a region string like 'chr1:100-200' into (chrom, start, end)."""
    match = REGION_REGEX.match(region_str.strip())
    if not match:
        raise ValueError(f"Malformed region string '{region_str}'. Expected 'chrom:start-end'.")
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def add_to_where_clause(query, where_clause):
    where_index = query.lower().find("where")
    if where_index == -1:
        # Check for ORDER BY, GROUP BY, LIMIT
        tail_match = re.search(r"\b(order\s+by|group\s+by|limit)\b", query, re.IGNORECASE)
        if tail_match:
            idx = tail_match.start()
            return f"{query[:idx]} WHERE {where_clause} {query[idx:]}"
        return f"{query} WHERE {where_clause}"
    prefix = query[:where_index]
    suffix = query[where_index + len("where") :]
    return f"{prefix} WHERE ({suffix.strip()}) AND ({where_clause})"


def add_region_to_query(query, region_str):
    """Add genomic coordinate filtering to a SQL query."""
    chrom, start, end = parse_region(region_str)
    where_clause = (
        f"chrom = '{chrom}' AND "
        f"((start BETWEEN {start} AND {end}) OR (end BETWEEN {start} AND {end}))"
    )
    return add_to_where_clause(query, where_clause)
