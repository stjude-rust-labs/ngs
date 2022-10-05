use std::{borrow::Borrow, cmp::Eq, collections::HashMap, hash::Hash, ops::Deref, sync::Mutex};

#[derive(Debug)]
pub struct MutexBackedHashMap<K, V>(Mutex<HashMap<K, V>>);

impl<K, V> MutexBackedHashMap<K, V>
where
    K: Eq + Hash + Clone,
    V: Clone,
{
    pub fn new(hashmap: HashMap<K, V>) -> Self {
        MutexBackedHashMap(Mutex::new(hashmap))
    }

    pub fn get<Q: ?Sized>(&self, k: &Q) -> Option<V>
    where
        K: Borrow<Q>,
        Q: Hash + Eq,
    {
        // SAFETY: MutexGuard dropped at the end of this scope.
        let hm = self.0.lock().unwrap();
        Some((*hm.get(k).unwrap()).clone())
    }

    pub fn insert(&mut self, k: K, v: V) -> Option<V> {
        // SAFETY: MutexGuard dropped at the end of this scope.
        let mut hm = self.0.lock().unwrap();
        hm.insert(k, v)
    }

    pub fn clone_nonatomic(&self) -> HashMap<K, V> {
        self.0.lock().unwrap().clone()
    }
}
