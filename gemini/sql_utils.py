import re


def get_select_cols_and_rest(query):
    """Separate the selected columns from the rest of the query.

    Returns:
        1. list of selected columns
        2. string of the rest of the query starting from FROM
    """
    from_loc = query.lower().find("from")
    if from_loc == -1:
        return [query.strip()], ""

    raw_select_clause = query[0:from_loc].rstrip()
    rest_of_query = query[from_loc:]

    select_pattern = re.compile(r"^\s*select\s+", re.IGNORECASE)
    raw_select_clause = select_pattern.sub("", raw_select_clause)

    selected_columns = [c.strip() for c in raw_select_clause.split(",") if c.strip()]
    return selected_columns, rest_of_query


def ensure_columns(query, cols):
    """Add columns to SELECT clause if they are not already present."""
    sel_cols, rest = get_select_cols_and_rest(query)
    if "*" in sel_cols:
        return query
    for c in cols:
        if c not in sel_cols:
            sel_cols.append(c)

    sel_string = ", ".join(sel_cols)
    return f"SELECT {sel_string} {rest}"
