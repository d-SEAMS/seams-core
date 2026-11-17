// Parse mermaid diagrams with the mermaid package named by MERMAID_PACKAGE.
// JSDOM_MODULE is the absolute path of jsdom's lib/api.js. Both are
// installed by docs/check_mermaid.py at the version in docs/source/conf.py.
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";

const mermaidRoot = process.env.MERMAID_PACKAGE;
const jsdomModule = process.env.JSDOM_MODULE;
if (!mermaidRoot || !jsdomModule) {
  console.error("MERMAID_PACKAGE and JSDOM_MODULE are required");
  process.exit(1);
}

const { JSDOM } = await import(pathToFileURL(jsdomModule).href);
const dom = new JSDOM("<!DOCTYPE html><html><body></body></html>");
globalThis.window = dom.window;
globalThis.document = dom.window.document;
globalThis.DOMParser = dom.window.DOMParser;
globalThis.HTMLElement = dom.window.HTMLElement;
globalThis.SVGElement = dom.window.SVGElement;
globalThis.Element = dom.window.Element;
globalThis.Node = dom.window.Node;

const expect = process.env.MERMAID_EXPECT;
const pkg = JSON.parse(readFileSync(`${mermaidRoot}/package.json`, "utf8"));
if (pkg.version !== expect) {
  console.error(`mermaid ${pkg.version} != conf.py ${expect}`);
  process.exit(1);
}

const mermaid = (
  await import(pathToFileURL(`${mermaidRoot}/dist/mermaid.core.mjs`).href)
).default;
await mermaid.initialize({ startOnLoad: false });

const blocks = JSON.parse(readFileSync(0, "utf8"));
let failed = 0;
for (const block of blocks) {
  try {
    await mermaid.parse(block.text);
  } catch (err) {
    failed += 1;
    const message = err && err.message ? err.message : String(err);
    console.error(`${block.path}\n${message}`);
  }
}
if (failed) {
  console.error(`${failed} mermaid diagram(s) failed`);
  process.exit(1);
}
console.log(`mermaid ${pkg.version} ok (${blocks.length} diagrams)`);
