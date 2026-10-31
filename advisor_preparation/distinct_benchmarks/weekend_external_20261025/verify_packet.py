from pathlib import Path
import hashlib, json, tarfile
p=Path(__file__).resolve().parent
m=json.loads((p/'PACKET_MANIFEST.json').read_text())
count=0
for f in m['auxiliary_files']:
 b=(p/f['path']).read_bytes();assert len(b)==f['bytes'];assert hashlib.sha256(b).hexdigest()==f['sha256']
for a in m['archives']:
 b=(p/a['path']).read_bytes();assert len(b)==a['bytes'];assert hashlib.sha256(b).hexdigest()==a['sha256']
 with tarfile.open(p/a['path'],'r:gz') as t:
  members=t.getmembers();assert len(members)==len(a['members']);assert [x.name for x in members]==sorted(x.name for x in members)
  assert len({x.name for x in members})==len(members)
  for x,d in zip(members,a['members']):
   assert x.isfile() and x.mode==0o644 and x.mtime==0
   assert not Path(x.name).is_absolute() and '..' not in Path(x.name).parts
   assert x.name==d['public_member'];v=t.extractfile(x).read();v.decode('utf-8')
   assert len(v)==d['public_bytes'] and hashlib.sha256(v).hexdigest()==d['public_sha256']
   if d['transformation']['type']=='BYTE_EXACT':assert d['original_sha256']==d['public_sha256'] and d['original_bytes']==d['public_bytes']
   count+=1
print(json.dumps({'status':'PASS','archives':len(m['archives']),'regular_UTF8_members':count,'scope':'TRANSPORT_INTEGRITY_ONLY_NO_SCIENTIFIC_REPLAY'},sort_keys=True))
