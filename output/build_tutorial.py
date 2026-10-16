"""Assemble the final offline tutorial; never run MD or modify research results.

Optional authoring command: .venv/Scripts/python.exe output/build_tutorial.py
Reading the delivered HTML requires no build command or companion resources.
"""
import base64
import html
import json
import re
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
SRC = OUT / 'tutorial-source'


def embedded(ident, value):
    text = json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
    return f'<script type="application/json" id="{ident}">' + text.replace('</', '<\\/') + '</script>'


def main():
    data = json.loads((OUT / 'tutorial-data.json').read_text(encoding='utf-8'))
    snippets = []
    for spec in json.loads((SRC / 'snippet-ranges.json').read_text(encoding='utf-8')):
        lines = (ROOT / spec['file']).read_text(encoding='utf-8').splitlines()
        snippets.append({'file': spec['file'], 'start': spec['start'],
                         'lines': lines[spec['start']-1:spec['start']-1+spec['count']]})
    figures = [{'caption': p.stem, 'data': 'data:image/png;base64,' + base64.b64encode(p.read_bytes()).decode()}
               for p in sorted((ROOT / 'results/analysis').glob('figure*.png'))]
    template = (SRC / 'template.html.in').read_text(encoding='utf-8')
    assert template.count('__TUTORIAL_PAYLOAD__') == 1
    ids = re.findall(r'\bid="([^"]+)"', template)
    assert len(ids) == len(set(ids)), 'Duplicate HTML IDs'
    assert len(re.findall(r'<section id=', template)) == 8
    assert not re.search(r'<(?:script|link)[^>]+(?:src|href)="https?://', template)
    script = (SRC / 'interactions.js').read_text(encoding='utf-8')
    payload = embedded('benchmark-data', data) + embedded('source-snippets', snippets) + embedded('figure-data', figures)
    page = template.replace('__TUTORIAL_PAYLOAD__', payload + '<script>' + script + '</script>')
    destination = OUT / 'ion-model-benchmark-tutorial-desktop.html'
    destination.write_text(page, encoding='utf-8')
    exercises = re.findall(r'<code data-practice-python="(\d+)">(.*?)</code>', template, re.S)
    assert len(exercises) == 6
    code = {number: html.unescape(value) for number, value in exercises}
    runner = '''"""Run a displayed answer: .venv/Scripts/python.exe output/tutorial-exercises.py 1
Choose 1..6 or all. Reads existing data; does not run production.
"""
from pathlib import Path
import os, sys
ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))
EXERCISES = ''' + repr(code) + '''
choice = sys.argv[1] if len(sys.argv) > 1 else "1"
if choice not in (*EXERCISES, "all"):
    raise SystemExit("Choose 1..6 or all")
for number in (EXERCISES if choice == "all" else [choice]):
    print("Exercise", number)
    exec(compile(EXERCISES[number], "exercise_" + number, "exec"), {"__name__": "__main__"})
'''
    (OUT / 'tutorial-exercises.py').write_text(runner, encoding='utf-8')
    print(f'Built {destination.name}: {destination.stat().st_size:,} bytes; 8 chapters; 6 Python exercises.')


if __name__ == '__main__':
    main()
