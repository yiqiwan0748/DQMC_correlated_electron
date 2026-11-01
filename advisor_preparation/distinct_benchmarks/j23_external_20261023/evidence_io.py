"""Streaming scientific UTF8 text archives, exact manifests and safe extraction."""
import codecs,gzip,hashlib,json,tarfile
from pathlib import Path,PurePosixPath
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def safe(name):
 p=PurePosixPath(name)
 if not name or p.is_absolute() or '..' in p.parts or '\\' in name or ':' in name:raise ValueError('Unsafe scientific member: '+name)
 return p
def text_info(p):
 p=Path(p)
 if p.is_symlink() or not p.is_file():raise ValueError('Only regular original text files are allowed')
 decode=codecs.getincrementaldecoder('utf-8')();h=hashlib.sha256();size=0
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):
   if b'\0' in b:raise ValueError('Runtime binary/NUL is not scientific text')
   decode.decode(b);h.update(b);size+=len(b)
 decode.decode(b'',final=True);return dict(sha256=h.hexdigest(),bytes=size)
def pack(path,entries):
 path=Path(path)
 if path.exists():raise ValueError('Fresh archive path required')
 members=[];fold=set();path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('xb') as output,gzip.GzipFile(fileobj=output,mode='wb',filename='',mtime=0) as compressed,tarfile.open(fileobj=compressed,mode='w|',format=tarfile.PAX_FORMAT) as tar:
  for name,p in sorted(entries.items()):
   safe(name)
   if name.casefold() in fold:raise ValueError('Cross-platform duplicate scientific member')
   fold.add(name.casefold());spec=text_info(p);info=tarfile.TarInfo(name);info.size=spec['bytes'];info.mode=0o644;info.uid=info.gid=0;info.mtime=0;info.uname=info.gname=''
   with Path(p).open('rb') as f:tar.addfile(info,f)
   if sha(p)!=spec['sha256']:raise ValueError('Original scientific input changed while archiving')
   members.append(dict(path=name,**spec))
 path.chmod(0o644);return dict(schema='scientific_text_archive_v1',archive=path.name,archive_sha256=sha(path),archive_bytes=path.stat().st_size,file_count=len(members),files=members,uncompressed_bytes=sum(x['bytes'] for x in members),source_bytes_changed=False,runtime_binaries_included=False)
def extract_verified(archive,manifest,destination):
 archive=Path(archive);destination=Path(destination)
 if destination.exists() or sha(archive)!=manifest['archive_sha256'] or archive.stat().st_size!=manifest['archive_bytes']:raise ValueError('Fresh destination and exact archive required')
 expected={r['path']:r for r in manifest['files']}
 if len(expected)!=manifest['file_count'] or len({n.casefold() for n in expected})!=len(expected):raise ValueError('Duplicate manifest paths')
 destination.mkdir(parents=True);seen=set()
 with tarfile.open(archive,'r:gz') as tar:
  for item in tar:
   name=item.name;relative=safe(name)
   if not item.isfile() or item.mode!=0o644 or name in seen or name not in expected or item.size!=expected[name]['bytes']:raise ValueError('Unexpected,linked,duplicate or wrong-sized member')
   target=destination.joinpath(*relative.parts)
   if not target.resolve().is_relative_to(destination.resolve()):raise ValueError('Extraction escape')
   target.parent.mkdir(parents=True,exist_ok=True);h=hashlib.sha256();decode=codecs.getincrementaldecoder('utf-8')()
   with tar.extractfile(item) as source,target.open('xb') as output:
    for b in iter(lambda:source.read(1048576),b''):
     if b'\0' in b:raise ValueError('Archive contains runtime binary/NUL')
     decode.decode(b);h.update(b);output.write(b)
   decode.decode(b'',final=True)
   if h.hexdigest()!=expected[name]['sha256']:raise ValueError('Scientific member hash differs')
   target.chmod(0o644);seen.add(name)
 if seen!=set(expected):raise ValueError('Incomplete scientific archive')
 return destination
