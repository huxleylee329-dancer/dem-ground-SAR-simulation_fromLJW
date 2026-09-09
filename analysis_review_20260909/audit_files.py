"""Read every original file without executing project code; build a review inventory."""
from pathlib import Path
import hashlib, json, re, struct, zlib, collections, sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
BINARY = {'.mat','.obj','.mexw64','.dll','.p','.jpg','.png','.tiff','.pdf'}
CLASSES = {1:'cell',2:'struct',3:'object',4:'char',5:'sparse',6:'double',7:'single',8:'int8',9:'uint8',10:'int16',11:'uint16',12:'int32',13:'uint32',14:'int64',15:'uint64'}

def element(data, off, endian):
    tag, size = struct.unpack_from(endian+'II',data,off)
    if tag >> 16:
        size=tag>>16; typ=tag&65535
        return typ,data[off+4:off+4+size],off+8
    return tag,data[off+8:off+8+size],off+8+size+(0 if tag==15 else (-size)%8)

def matmeta(data):
    if b'MATLAB 7.3' in data[:128]: return {'format':'MATLAB 7.3/HDF5','status':'header only'}
    if not data.startswith(b'MATLAB 5.0'): return {'format':'other','status':'header only'}
    endian='<' if data[126:128]==b'IM' else '>'
    variables=[]
    def walk(buf, start):
        while start+8<=len(buf):
            typ,payload,nxt=element(buf,start,endian)
            if nxt<=start or nxt>len(buf)+7: raise ValueError('Invalid MAT element')
            start=nxt
            if typ==15: walk(zlib.decompress(payload),0)
            if typ==14:
                _,flags,q=element(payload,0,endian)
                _,dims,q=element(payload,q,endian)
                _,name,q=element(payload,q,endian)
                f=struct.unpack_from(endian+'I',flags)[0]
                shape=struct.unpack(endian+'i'*(len(dims)//4),dims)
                v={'name':name.decode('utf-8','replace'),'shape':shape,'class':CLASSES.get(f&255,str(f&255)),'complex':bool(f&0x800)}
                if f&255 in range(6,16) and len(payload)>q:
                    t,d,_=element(payload,q,endian)
                    formats={1:'b',2:'B',3:'h',4:'H',5:'i',6:'I',7:'f',9:'d',12:'q',13:'Q'}
                    if t in formats and 0<len(d)<=128:
                        fmt=formats[t]; v['small_values']=struct.unpack(endian+fmt*(len(d)//struct.calcsize(fmt)),d)
                variables.append(v)
    walk(data,128)
    return {'format':'MATLAB 5 compressed/uncompressed','variables':variables}

records=[]; texts={}; failures=[]
for path in sorted(ROOT.rglob('*')):
    if not path.is_file() or OUT in path.parents: continue
    rel=path.relative_to(ROOT).as_posix()
    try:
        data=path.read_bytes()
        item={'path':rel,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
        if path.suffix.lower() not in BINARY:
            for enc in ['utf-8-sig','gb18030','cp1252']:
                try: s=data.decode(enc); break
                except UnicodeDecodeError: pass
            else: s=data.decode('utf-8','replace');enc='utf-8 replacement'
            texts[rel]=s
            lines=s.splitlines()
            item.update(kind='text',encoding=enc,lines=len(lines))
            item['declarations']=[l.strip() for l in lines if re.match(r'^\s*(function|classdef)\b',l)][:12]
            item['comments']=[l.strip()[:180] for l in lines if l.strip().startswith(('%','#','//'))][:5]
            item['sections']=[l.strip()[:180] for l in lines if l.strip().startswith('%%')][:30]
            item['io']=[{'line':i+1,'text':l.strip()[:250]} for i,l in enumerate(lines) if not l.lstrip().startswith('%') and re.search(r'\b(load|save|csvread|fopen|imread|addpath|readbinary|writebinary)\b',l)]
        else:
            item.update(kind='binary',signature=data[:24].hex())
            if path.suffix.lower()=='.mat':
                try: item['mat']=matmeta(data)
                except Exception as e: item['mat']={'error':str(e)}
        records.append(item)
    except Exception as e: failures.append({'path':rel,'error':str(e)})
(OUT/'inventory.json').write_text(json.dumps({'files':records,'failures':failures},ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'source_texts.json').write_text(json.dumps(texts,ensure_ascii=False),encoding='utf-8')
catalog=['# 逐文件读取清单','', '源码和文本完整读取；二进制完整读取字节并计算 SHA-256，MAT 尽可能解析顶层变量。二进制校验不等于反编译或逐像素审阅。','']
for item in records:
    details=f"{item['bytes']} bytes"
    if item['kind']=='text':
        details+=f"; {item['lines']} lines; {item['encoding']}"
        details+='; '+ ' / '.join(item['declarations'][:1]+item['comments'][:2])
    elif 'mat' in item:
        details+='; '+json.dumps(item['mat'],ensure_ascii=False)
    else: details+='; 二进制字节读取，未反编译/渲染'
    catalog.append('- `'+item['path']+'` — '+details.replace('\n',' '))
(OUT/'逐文件读取清单.md').write_text('\n'.join(catalog),encoding='utf-8')
print(json.dumps({'files':len(records),'bytes':sum(i['bytes'] for i in records),'text_files':len(texts),'text_lines':sum(i.get('lines',0) for i in records),'binary_files':sum(i['kind']=='binary' for i in records),'failures':failures,'mat_formats':collections.Counter(i['mat'].get('format','error') for i in records if 'mat' in i)},ensure_ascii=False))
