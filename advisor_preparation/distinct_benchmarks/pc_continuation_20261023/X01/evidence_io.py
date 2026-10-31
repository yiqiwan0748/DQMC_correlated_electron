"""Portable, byte-exact scientific text archives; no solver execution."""
from pathlib import Path, PurePosixPath
import gzip,hashlib,io,json,tarfile

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def safe(name):
    p=PurePosixPath(name)
    if not name or p.is_absolute() or '..' in p.parts or '\\' in name or ':' in name:
        raise ValueError('Unsafe archive path: '+name)
    return p

def pack(path, entries):
    """entries maps explicit safe relative names to immutable bytes."""
    path=Path(path);seen=set();members=[]
    if path.exists():raise ValueError('Use a new archive')
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as output,gzip.GzipFile(fileobj=output,mode='wb',filename='',mtime=0) as compressed,tarfile.open(fileobj=compressed,mode='w|',format=tarfile.PAX_FORMAT) as tar:
        for name,value in sorted(entries.items()):
            safe(name)
            if name.casefold() in seen:raise ValueError('Cross-platform duplicate member')
            seen.add(name.casefold());value.decode('utf-8')
            if b'\0' in value:raise ValueError('Runtime binary is not scientific text')
            info=tarfile.TarInfo(name);info.size=len(value);info.mode=0o644;info.uid=info.gid=0;info.mtime=0;info.uname=info.gname=''
            tar.addfile(info,io.BytesIO(value));members.append({'path':name,'bytes':len(value),'sha256':hashlib.sha256(value).hexdigest()})
    path.chmod(0o644)
    return {'schema':'scientific_text_archive_v1','archive':path.name,'archive_sha256':sha(path),'archive_bytes':path.stat().st_size,'files':members,'file_count':len(members),'uncompressed_bytes':sum(r['bytes'] for r in members),'source_bytes_changed':False,'runtime_binaries_included':False}

def extract_verified(archive, manifest, destination):
    archive=Path(archive);destination=Path(destination)
    if destination.exists():raise ValueError('Use a fresh extraction directory')
    if sha(archive)!=manifest['archive_sha256'] or archive.stat().st_size!=manifest['archive_bytes']:
        raise ValueError('Archive digest differs')
    expected={r['path']:r for r in manifest['files']}
    if len(expected)!=manifest['file_count'] or len({n.casefold() for n in expected})!=len(expected):raise ValueError('Duplicate manifest paths')
    destination.mkdir(parents=True);seen=set()
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar:
            relative=safe(member.name)
            if not member.isfile() or member.name in seen or member.name not in expected:raise ValueError('Unexpected/link/duplicate member')
            row=expected[member.name]
            if member.size!=row['bytes']:raise ValueError('Unexpected member size')
            value=tar.extractfile(member).read()
            if hashlib.sha256(value).hexdigest()!=row['sha256']:raise ValueError('Member digest differs')
            target=destination.joinpath(*relative.parts)
            if not target.resolve().is_relative_to(destination.resolve()):raise ValueError('Extraction escape')
            target.parent.mkdir(parents=True,exist_ok=True)
            with target.open('xb') as out:out.write(value)
            target.chmod(0o644);seen.add(member.name)
    if seen!=set(expected):raise ValueError('Incomplete archive')
    return destination
