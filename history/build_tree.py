"""Rebuild the public, standalone tree from the included JSON and template."""
from pathlib import Path
import json

root = Path(__file__).resolve().parent
data = json.loads((root / 'research_tree_data.json').read_text(encoding='utf-8'))
template = (root / 'tree_template.html').read_text(encoding='utf-8')
assert template.count('__RESEARCH_DATA__') == 1
payload = json.dumps(data, ensure_ascii=False).replace('</', '<\\/')
output = root / 'research-evolution-tree.html'
output.write_text(template.replace('__RESEARCH_DATA__', payload), encoding='utf-8')
print(f"Built {output.name}: {len(data['nodes'])} records, {len(data['stages'])} stages")
