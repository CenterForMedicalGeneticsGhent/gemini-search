import sys

from gemini import gemini_query


def query_gemini(db: str, query: str, gt_query: str | None = None):
    gq = gemini_query.GeminiQuery(
        db,
        out_format=gemini_query.JSONRowFormat(None),
        include_gt_cols=False,
    )
    if gt_query:
        gq.run(query, gt_query)
    else:
        gq.run(query)
    return gq


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python test_full_db_query.py <db> <query> [<gt_query>]")
        sys.exit(1)
    db = sys.argv[1]
    query = sys.argv[2]
    gt_query = sys.argv[3] if len(sys.argv) > 3 else None
    result = query_gemini(db, query, gt_query)
    rows = [{k: row[k] for k in row} for row in result]
    print(f"Retrieved {len(rows)} rows.")
