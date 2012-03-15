"""
Import the content of one or more files to the database.

This can be used to import even big Magnet and RF files into the database. 
The files still need to be linked with the accelerator that uses them. 
"""
import os
import sys

from util.database import DatabaseConnection

def main():
    chunk_size = 1024*512 # 512kb
    
    if len(sys.argv) < 2:
        sys.exit('Syntax: file_import.py <path, path, ...>')
    
    files = sys.argv[1:]
    for path in files:
        if not os.path.isfile(path):
            print 'File not found: %s, omitting!' % path
            files.remove(path)
    
    if len(files) < 1:
        sys.exit('... No files to push, good bye!')
    
    conn = DatabaseConnection()
    for path in files:
        print 'Pushing %s to the database...' % path
            
        file_id = conn.field("SELECT id FROM file WHERE filename = '%s'" % os.path.basename(path))
        if file_id is None:
            print '... no such file in the database, creating new entry'
            file_id = conn.insert('file', ['content', 'filename'], [('', os.path.basename(path))])
        else:
            print '... found file in the database, overwriting content'
            conn.execute("UPDATE file SET content = '' WHERE id = %s" % file_id)
        
        query = """ UPDATE file SET content = CONCAT(content, '%(content)s') WHERE id = %(id)s """    
        with open(path, 'r') as f:
            print '... opened %s for reading' % path
            sys.stdout.write('... writing contents to the database: ')
            while True:
                chunk = f.read(chunk_size)
                if len(chunk) == 0:
                    break
                conn.execute(query % {'id': file_id, 'content': conn.escape_string(chunk)})
                sys.stdout.write('.')
                sys.stdout.flush()
        
        print '... done'                

if __name__ == '__main__':
    main()
