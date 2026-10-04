// generate.js — 论文 Word 文档渲染器（预印本样式：标题块 + 摘要 + 目录 + 正文）
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  HeadingLevel, AlignmentType, WidthType, BorderStyle, ShadingType,
  TableOfContents, PageBreak, Footer, Header, PageNumber, NumberFormat, SectionType,
} = require("docx");
const fs = require("fs");
const { C } = require("./content.js");
const { C2 } = require("./content2.js");

const FONT = { eastAsia: "SimSun", ascii: "Times New Roman" };
const FHEI = { eastAsia: "SimHei", ascii: "Times New Roman" };
const FTIME = { ascii: "Times New Roman", eastAsia: "SimSun" };
const BLACK = "000000";
const NB = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const noBorders = { top: NB, bottom: NB, left: NB, right: NB };

// 行内解析：_{..} 下标、^{..} 上标
function parseRuns(text, base) {
  const runs = [];
  const re = /([_^])\{([^}]*)\}/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) runs.push(new TextRun({ ...base, text: text.slice(last, m.index) }));
    if (m[1] === "_") runs.push(new TextRun({ ...base, text: m[2], subScript: true }));
    else runs.push(new TextRun({ ...base, text: m[2], superScript: true }));
    last = re.lastIndex;
  }
  if (last < text.length) runs.push(new TextRun({ ...base, text: text.slice(last) }));
  if (runs.length === 0) runs.push(new TextRun({ ...base, text }));
  return runs;
}

const bodyBase = { size: 24, color: BLACK, font: FONT };
const fmBase = { size: 24, color: BLACK, font: FTIME };

function bodyPara(text, opts = {}) {
  return new Paragraph({
    alignment: AlignmentType.JUSTIFIED,
    spacing: { line: 312, before: 40, after: 40 },
    indent: opts.noIndent ? undefined : { firstLine: 480 },
    children: parseRuns(text, bodyBase),
  });
}

function formulaPara(text) {
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { line: 312, before: 100, after: 100 },
    children: parseRuns(text, { ...fmBase, italics: false }),
  });
}

function thmPara(text) {
  return new Paragraph({
    alignment: AlignmentType.LEFT,
    spacing: { line: 312, before: 160, after: 60 },
    keepNext: true,
    children: [new TextRun({ text, bold: true, size: 24, color: BLACK, font: FHEI })],
  });
}

function notePara(text) {
  return new Paragraph({
    alignment: AlignmentType.JUSTIFIED,
    spacing: { line: 312, before: 40, after: 40 },
    children: parseRuns(text, { size: 21, color: "555555", font: FONT, italics: true }),
  });
}

function heading(text, level, opts = {}) {
  const conf = {
    1: { size: 32, bold: true, centered: true, before: 300, after: 160 },
    2: { size: 28, bold: true, centered: false, before: 240, after: 120 },
    3: { size: 24, bold: true, centered: false, before: 200, after: 100 },
  }[level];
  return new Paragraph({
    heading: level === 1 ? HeadingLevel.HEADING_1 : (level === 2 ? HeadingLevel.HEADING_2 : HeadingLevel.HEADING_3),
    alignment: conf.centered ? AlignmentType.CENTER : AlignmentType.LEFT,
    spacing: { line: 312, before: conf.before, after: conf.after },
    keepNext: true,
    children: [new TextRun({ text, size: conf.size, bold: conf.bold, color: BLACK, font: FHEI })],
  });
}

const cellMargins = { top: 60, bottom: 60, left: 100, right: 100 };
function buildTable(spec) {
  const out = [];
  if (spec.caption) {
    out.push(new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { before: 160, after: 80, line: 312 },
      keepNext: true,
      children: [new TextRun({ text: spec.caption, bold: true, size: 21, color: BLACK, font: FHEI })],
    }));
  }
  const headerRow = new TableRow({
    tableHeader: true,
    cantSplit: true,
    children: spec.headers.map((h, i) => new TableCell({
      width: { size: spec.widths[i], type: WidthType.PERCENTAGE },
      shading: { type: ShadingType.CLEAR, fill: "F5F7FA", color: "auto" },
      margins: cellMargins,
      children: [new Paragraph({
        alignment: AlignmentType.CENTER, spacing: { line: 280 },
        children: [new TextRun({ text: h, bold: true, size: 21, color: BLACK, font: FHEI })],
      })],
    })),
  });
  const rows = spec.rows.map(r => new TableRow({
    cantSplit: true,
    children: r.map((c, i) => new TableCell({
      width: { size: spec.widths[i], type: WidthType.PERCENTAGE },
      margins: cellMargins,
      children: [new Paragraph({
        alignment: i === 0 ? AlignmentType.LEFT : AlignmentType.LEFT,
        spacing: { line: 280 },
        children: parseRuns(c, { size: 21, color: BLACK, font: FONT }),
      })],
    })),
  }));
  out.push(new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    borders: {
      top: { style: BorderStyle.SINGLE, size: 6, color: BLACK },
      bottom: { style: BorderStyle.SINGLE, size: 6, color: BLACK },
      left: NB, right: NB,
      insideHorizontal: { style: BorderStyle.SINGLE, size: 2, color: "999999" },
      insideVertical: NB,
    },
    rows: [headerRow, ...rows],
  }));
  out.push(new Paragraph({ spacing: { after: 120 }, children: [] }));
  return out;
}

function render(blocks) {
  const out = [];
  for (const b of blocks) {
    if (b.h1 !== undefined) out.push(heading(b.h1, 1));
    else if (b.h2 !== undefined) out.push(heading(b.h2, 2));
    else if (b.h3 !== undefined) out.push(heading(b.h3, 3));
    else if (b.thm !== undefined) out.push(thmPara(b.thm));
    else if (b.pf !== undefined) out.push(bodyPara(b.pf));
    else if (b.fm !== undefined) out.push(formulaPara(b.fm));
    else if (b.p !== undefined) out.push(bodyPara(b.p));
    else if (b.pn !== undefined) out.push(bodyPara(b.pn, { noIndent: true }));
    else if (b.note !== undefined) out.push(notePara(b.note));
    else if (b.tbl !== undefined) out.push(...buildTable(b.tbl));
  }
  return out;
}

/* ---------- 第一节：标题块 + 摘要 + 目录 ---------- */
const front = [];
front.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { before: 600, after: 120, line: 340 },
  children: [new TextRun({ text: "偏序集 1/3–2/3 猜想的一类 width-3 情形", bold: true, size: 36, color: BLACK, font: FHEI })],
}));
front.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { after: 120, line: 340 },
  children: [new TextRun({ text: "完整归约、可证工具箱与守恒纲领", bold: true, size: 30, color: BLACK, font: FHEI })],
}));
front.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { after: 60, line: 280 },
  children: [new TextRun({ text: "工作论文（研究札记）", size: 21, color: "555555", font: FONT })],
}));
front.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { after: 300, line: 280 },
  children: [new TextRun({ text: "整理自 2026 年 9 月至 10 月的攻击会话 · 档案可复跑", size: 21, color: "555555", font: FONT })],
}));
// 摘要（标题块后直接排）
front.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { before: 200, after: 120, line: 312 }, keepNext: true,
  children: [new TextRun({ text: "摘  要", bold: true, size: 28, color: BLACK, font: FHEI })],
}));
for (const b of C.filter(x => x.h1 === "摘要")) {
  if (b.p !== undefined) front.push(bodyPara(b.p));
}
front.push(new Paragraph({
  alignment: AlignmentType.LEFT, spacing: { before: 300, after: 120, line: 312 }, keepNext: true,
  children: [new TextRun({ text: "目  录", bold: true, size: 28, color: BLACK, font: FHEI })],
}));
front.push(new TableOfContents("目录", { hyperlink: true, headingStyleRange: "1-2" }));
front.push(new Paragraph({
  spacing: { before: 160 },
  children: [new TextRun({
    text: "注：目录由域代码生成，编辑文档后请右键目录并选择更新域以刷新页码。",
    italics: true, size: 18, color: "888888", font: FONT,
  })],
}));

/* ---------- 第二节：正文 ---------- */
const body = render(C.filter(x => x.h1 !== "摘要").concat(C2));

/* ---------- 文档 ---------- */
const doc = new Document({
  creator: "poset13-attack",
  title: "偏序集 1/3–2/3 猜想的一类 width-3 情形",
  styles: {
    default: {
      document: { run: { size: 24, color: BLACK, font: FONT } },
    },
  },
  features: { updateFields: true },
  sections: [
    {
      properties: {
        page: {
          size: { width: 11906, height: 16838 },
          margin: { top: 1440, bottom: 1440, left: 1701, right: 1417, header: 850, footer: 992 },
        },
      },
      headers: {
        default: new Header({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [new TextRun({ text: "偏序集 1/3–2/3 猜想的一类 width-3 情形 · 工作论文", size: 18, color: "333333", font: FONT })],
          })],
        }),
      },
      footers: {
        default: new Footer({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [new TextRun({ children: [PageNumber.CURRENT], size: 21, color: BLACK, font: { ascii: "Times New Roman" } })],
          })],
        }),
      },
      children: front,
    },
    {
      properties: {
        type: SectionType.NEXT_PAGE,
        page: {
          size: { width: 11906, height: 16838 },
          margin: { top: 1440, bottom: 1440, left: 1701, right: 1417, header: 850, footer: 992 },
          pageNumbers: { start: 1, formatType: NumberFormat.DECIMAL },
        },
      },
      headers: {
        default: new Header({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [new TextRun({ text: "偏序集 1/3–2/3 猜想的一类 width-3 情形 · 工作论文", size: 18, color: "333333", font: FONT })],
          })],
        }),
      },
      footers: {
        default: new Footer({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [new TextRun({ children: [PageNumber.CURRENT], size: 21, color: BLACK, font: { ascii: "Times New Roman" } })],
          })],
        }),
      },
      children: body,
    },
  ],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync("paper.docx", buf);
  console.log("OK paper.docx", buf.length, "bytes");
});
